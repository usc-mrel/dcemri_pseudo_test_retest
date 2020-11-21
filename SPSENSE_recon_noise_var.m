function [img,opt] = SPSENSE_recon_noise_var(k,opt)
%   CS-SENSE reconstruction
%   
%   Author: RML
%   Date: 08/2011
%   
%   Usage: [img,opt] = SPSENSE_recon(kU,opt)
%   
%   Input:
%   kU: Undersampled k-space data of size RO x PE x NS x NT x NR
%       NB: must have missing lines filled with zeros
%           should not be fftshifed (i.e., power should be at corners)
%   opt: Specifies reconstruction options. See SPSENSE_optset.m
%   
%   Output:
%   img: reconstructed image
%   opt: modified reconstruction options (contains PI calibration, k-space
%       weights, reconstruction info, etc...)


%   Start the recon timer
ReconTime = tic;
opt.RI.ReconTime = toc(ReconTime);

%   Check inputs
if nargin < 2
    error('SPSENSE_recon requires at least two inputs');
end

%   Get and save data size
np = opt.size(1); nv = opt.size(2); ns = opt.size(3);
nt = opt.size(4); nr = opt.size(5);
FTnorm = prod(opt.size(opt.FTdim));
%FTnorm = 1;

%   Check to ensure requested parameters are consistent with input data
%   Disable GRAPPA (not used in SENSE)
%   Disable sparcity transforms along singleton dimensions
param_check;

%   Check input data types to ensure correct class
class_check;

%   Define optimization parameters for line search
opt2 = optimset('fminbnd');
opt2.Display = 'off';
opt2.TolX = 1e-3;
opt2.MaxFunEvals = 40;
opt2.MaxIter = 40;
tmax = 4;

%   Estimate initial image
comp_img0;

%   Estimate image phase and apply
comp_Ph;
img = opt.iPH(img,opt);

%   Estimate wavelet transform order (if not specified)
comp_worder;

%   Initialize computation variables
knew = 0; kd   = 0;
WC   = 0; WCd  = 0;
TV1  = 0; TV1d = 0;
TV2  = 0; TV2d = 0;
VS   = 0; VSd  = 0;
RI   = 0; RId  = 0;
CD   = 0; CDd  = 0;
if ~isfield(opt,'WCwt') || isempty(opt.WCwt)
    opt.WCwt = 0;
end
if ~isfield(opt,'TV1wt') || isempty(opt.TV1wt)
    opt.TV1wt = 0;
end
if ~isfield(opt,'TV2wt') || isempty(opt.TV2wt)
    opt.TV2wt = 0;
end
if ~isfield(opt,'VSwt') || isempty(opt.VSwt)
    opt.VSwt = 0;
end
if ~isfield(opt,'RIwt') || isempty(opt.RIwt)
    opt.RIwt = 0;
end
if ~isfield(opt,'CDwt') || isempty(opt.CDwt)
    opt.CDwt = 0;
end

%   Initialize l1-weights
comp_l1weight;

%   Compute initial gradient
G = 0;
CompGrad;
delta = -G;
nrmG = G(:)'*G(:);

%   Compute initial objective function
SPIRiTpreobj;
[opt.RI.Obj, opt.RI.DC,  opt.RI.SP1, opt.RI.SP2, ...
 opt.RI.SP3, opt.RI.SP4, opt.RI.SP5, opt.RI.SP6] = SPIRiT_objective(0);
opt.RI.nrmG = nrmG;

%   Plot, if requested
opt.RI.Iter = 0;
if opt.plot
    close all;
    opt = opt.plotf(img,opt);
end
if opt.verbose == 2
    fprintf('Iter.: %3.0u Obj: %7.5f\n',[0 opt.RI.Obj]);
end

%   Iterate through conjugate gradient steps
Dobj = ones(1,opt.Dobjn);
% imgiter = [];
while opt.RI.Iter < opt.MinIter || (opt.RI.Iter < opt.MaxIter && median(Dobj) > opt.Dobj)
        
    %   Store subset of image
%     imgiter = cat(3,imgiter,squeeze(img(:,:,round(ns/2),round(nt/2))));
    
    %   Save previous objective function
    %   Compute pre-objective values to accelerate the objective function
    if opt.RI.Iter == 0
        f0 = opt.RI.Obj;
    else
        f0 = f1;
        SPIRiTpreobj;
    end
    
    %   Line search for minimum objective
    [t,f1] = fminbnd(@SPIRiT_objective,eps,tmax,opt2);
    if t > 0.5*tmax
        tmax = 5*tmax;
    elseif t < 0.01*tmax
        tmax = 0.5*tmax;
    end
    opt2.TolX = 1e-1*t;
    
    if t < 0 && opt.verbose
        disp('t is negative');
        opt2.TolX = 1e-5*tmax;
    end
%     if f1 > f0 && opt.verbose
%         disp('Wrong direction!');
%     end
    
    %   Update image with optimal change
%     img = mult_add_MEX(1,img,t,delta);
    img = img+t*delta;
    
    %   Update iteration counter
    opt.RI.Iter = opt.RI.Iter+1;
    
    %   Track changes in the objective function
    Dobj = circshift(Dobj,[0 -1]);
    Dobj(end) = abs(100*(f0-f1)./f1);
    
    %   Record penalties
    if nargout > 1
        [f1,dct,sp1,sp2,sp3,sp4,sp5,sp6] = SPIRiT_objective(t);
        opt.RI.Obj = [opt.RI.Obj f1];
        opt.RI.DC  = [opt.RI.DC  dct];
        opt.RI.SP1 = [opt.RI.SP1 sp1];
        opt.RI.SP2 = [opt.RI.SP2 sp2];
        opt.RI.SP3 = [opt.RI.SP3 sp3];
        opt.RI.SP4 = [opt.RI.SP4 sp4];
        opt.RI.SP5 = [opt.RI.SP5 sp5];
        opt.RI.SP6 = [opt.RI.SP6 sp6];
    end

    if opt.verbose == 2
        fprintf('It.: %3.0u; Obj: %7.5f; D_Obj: %7.5f; |G|: %7.5f\n',...
            [opt.RI.Iter f1 100*(f0-f1)/f1],nrmG);
    end
    
%     if nr > 1
%         figure(1210);images(abs(tile_2dwt(opt.fSP1(sos(img),opt.wname,opt.worder))));drawnow;
%     else
%         figure(1210);images(abs(tile_2dwt(opt.fSP1(img,opt.wname,opt.worder))));
%         if ~isempty(WCwt)
%             figure(1211);images(abs(tile_2dwt(WCwt)));
%         end
%         figure(1);drawnow;
%     end
    
    %   Update periodic items
    if opt.plot && opt.RI.Iter == 1
        opt = opt.plotf(img,opt);
    end
    if ~mod(opt.RI.Iter,opt.update)
        if opt.plot
            opt = opt.plotf(img,opt);
        end
        PeriodicUpdate;
    end
    
    %   Free memory prior to gradient calculation
    knew = 0; kd   = 0;
    WC   = 0; WCd  = 0;
    TV1  = 0; TV1d = 0;
    TV2  = 0; TV2d = 0;
    VS   = 0; VSd  = 0;
    RI   = 0; RId  = 0;
    CD   = 0; CDd  = 0;
    
    %   Update gradient and search direction (conjugated gradient)
%     CompGrad;
%     G_new = G(:)'*G(:);
%     gamma = G_new/(G_old+eps);
%     delta = mult_add_MEX(gamma,delta,-1,G); % same as "delta = gamma*delta - G;"
% %     delta = gamma*delta - G;
%     G_old = G_new;
%     G = 0;
    
    G_old = G;
    CompGrad;
    nrmG = G(:)'*G(:);
    opt.RI.nrmG = [opt.RI.nrmG nrmG];
    num = (G(:)'*(G(:)-G_old(:)));
    if opt.cg_mtd == 0
        dem = -(G_old(:)'*(G(:)-G_old(:)));
    else
        dem = (G_old(:)'*G_old(:));
    end
    gamma = real(num/dem);
    if gamma < 0
        grad_restart;
        opt.sqrt_smooth = 0.1*opt.sqrt_smooth;
        opt.cg_mtd = ~opt.cg_mtd;
        delta = -G;
    else
%         delta = mult_add_MEX(gamma,delta,-1,G); % same as "delta = gamma*delta - G;"
        delta = gamma*delta - G;
    end
    
    %   Save recon time
    opt.RI.ReconTime = toc(ReconTime);
    
end

%   Save intermediate images
% save('ImgIter.mat','imgiter');

%   Reapply phase
img = opt.fPH(img,opt);

%   Save some recon information
opt.RI.ReconTime = toc(ReconTime);

%   Display some recon info
if opt.verbose
    fprintf('Final constraints: l(1) = %6.5f; l(2) = %6.5f; l(3) = %6.5f; l(4) = %6.5f; l(5) = %6.5f; l(6) = %6.5f; l(7) = %6.5f\n',...
        opt.lambda);
end

%   Get recon time
if opt.verbose
    fprintf('Recon time: %6.2f minutes\n\n',opt.RI.ReconTime/60);
end





%   ------------------------------------------------------------------  %
%   Begin sub-functions (nested for simplicity and memory efficienty)   %
%   -----------------------------------------------------------------   %

%   Compute gradient
    function CompGrad
        
        %   Compute current k-space
        knew = opt.fFTU(opt.fPI(opt.fPH(img,opt),opt),opt);
        
        %   Compute data consistency gradient
        kdiff = (knew(:) - k(:));
        %kdiff = opt.kW(:).*(knew - k(:));
        G = 2*opt.iPH(opt.iPI(opt.iFTU(kdiff,opt),opt),opt);
        clear kdiff;
        
        %   Compute wavelet gradient
        if opt.lambda(2)
            WC = opt.fSP1(img,opt.wname,opt.worder);
            G = G + opt.lambda(2)*opt.iSP1(opt.fWT(opt.WCwt,absgrad(WC,opt.sqrt_smooth)),opt.wname,opt.worder);
        end
        
        %   Compute TV gradient
        if opt.lambda(3)
            TV1 = opt.fSP2(img,opt.FDorder1);
            G = G + opt.lambda(3)*opt.iSP2(opt.fWT(opt.TV1wt,absgrad(TV1,opt.sqrt_smooth)),opt.FDorder1);
        end
        
        %   Compute spatial k-space filter
        if opt.lambda(4)
            TV2 = opt.fSP3(img,opt);
            G = G + opt.lambda(4)*opt.iSP3(opt.fWT(opt.TV2wt,absgrad(TV2,opt.sqrt_smooth)),opt);
        end
        
        %   Compute view-sharing gradient
        if opt.lambda(5)
            VS = opt.fVS(img,opt);
            G = G + opt.lambda(5)*opt.iVS(opt.fWT(opt.VSwt,absgrad(VS,opt.sqrt_smooth)),opt);
        end
        
        %   Compute reference image gradient
        if opt.lambda(6)
            RI = opt.fSP4(img,opt);
            G = G + opt.lambda(6)*opt.fWT(opt.RIwt,absgrad(RI,opt.sqrt_smooth));
        end
        
        %   Compute composite dynamic gradient
        if opt.lambda(7)
            CD = opt.fCD(img,opt);
%             G = G + opt.lambda(7)*opt.iCD(opt.fWT(opt.CDwt,absgrad(CD,opt.sqrt_smooth)),opt);
            G = G + (2*opt.lambda(7))*opt.iCD(CD,opt);
        end
        
    end



%   Compute pre-objective function values
    function SPIRiTpreobj
        
        %   Compute fourier domain
        kd = opt.fFTU(opt.fPI(opt.fPH(delta,opt),opt),opt);
        
        %   Compute wavelet domain
        if opt.lambda(2)
            WCd = opt.fSP1(delta,opt.wname,opt.worder);
            WCd = opt.fWT(opt.WCwt,WCd);
            WC  = opt.fWT(opt.WCwt,WC);
        end
        
        %   Compute 1st total variation domain
        if opt.lambda(3)
            TV1d = opt.fSP2(delta,opt.FDorder1);
            TV1d = opt.fWT(opt.TV1wt,TV1d);
            TV1  = opt.fWT(opt.TV1wt,TV1);
        end
        
        %   Compute 2nd total variation domain
        if opt.lambda(4)
            TV2d = opt.fSP3(delta,opt);
            TV2d = opt.fWT(opt.TV2wt,TV2d);
            TV2  = opt.fWT(opt.TV2wt,TV2);
        end
        
        %   Compute view sharing domain
        if opt.lambda(5)
            VSd = opt.fVS(delta,opt);
            VSd = opt.fWT(opt.VSwt,VSd);
            VS  = opt.fWT(opt.VSwt,VS);
        end
        
        %   Compute refernce image domain
        if opt.lambda(6)
            RId = delta;
            RId = opt.fWT(opt.RIwt,RId);
            RI  = opt.fWT(opt.RIwt,RI);
        end
        
        %   Compute composite dynamic domain
        if opt.lambda(7)
            CDd = opt.fCD(delta,opt);
            CDd = opt.fWT(opt.RIwt,CDd);
            CD  = opt.fWT(opt.RIwt,CD);
        end
        
    end



    %   Objective function
    function [Objtve,DCf,SP1f,SP2f,SP3f,VSf,RIf,CDf] = SPIRiT_objective(t)
        
        %   Compute data consistency
        DCf = l2_norm_k_MEX(knew,kd,k,t) / FTnorm;
%         DCf = sum(abs(knew(:)+t*kd(:) - k(:)).^2) / FTnorm;
%         DCf = sum((opt.kW(:).*abs(knew+t*kd - k(:))).^2) / FTnorm;
        
        %   Compute wavelet norm
        if opt.lambda(2)
            SP1f = opt.lambda(2) * l1norm(WC,WCd,t,opt.sqrt_smooth);
        else
            SP1f = 0;
        end
        
        %   Compute 1st TV norm
        if opt.lambda(3)
            SP2f = opt.lambda(3) * l1norm(TV1,TV1d,t,opt.sqrt_smooth);
        else
            SP2f = 0;
        end
        
        %   Compute 2nd TV norm
        if opt.lambda(4)
            SP3f = opt.lambda(4) * l1norm(TV2,TV2d,t,opt.sqrt_smooth);
        else
            SP3f = 0;
        end
        
        %   Compute view-sharing norm
        if opt.lambda(5)
            VSf = opt.lambda(5) * l1norm(VS,VSd,t,opt.sqrt_smooth);
        else
            VSf = 0;
        end
        
        %   Compute reference image norm
        if opt.lambda(6)
            RIf = opt.lambda(6) * l1norm(RI,RId,t,opt.sqrt_smooth);
        else
            RIf = 0;
        end
        
        %   Compute composite dynamic norm
        if opt.lambda(7)
%             CDf = opt.lambda(7) * l1norm(CD,CDd,t,opt.sqrt_smooth);
            CDf = opt.lambda(7) * sum(abs(CD(:)+t*CDd(:)).^2);
        else
            CDf = 0;
        end
        
        %   Compute net objective
        Objtve = DCf + SP1f + SP2f + SP3f + VSf + RIf + CDf;
    end



%   Periodic update function
%   This function adjusts opt.lambda to approach the target SNR
%   It also resets the conjugate gradient direction
    function PeriodicUpdate
        
        %   Flag to force restart conjugate gradient search direction
        grad_restart_flag = 0;
        
        %   Compute relative contributions to net objective function
        [Objtve,DCf,SP1f,SP2f,SP3f,VSf,RIf,CDf] = SPIRiT_objective(t);
        if opt.verbose == 2
            fprintf('-------\nObj: %6.2f; DC: %6.3f',Objtve,DCf);
            for i = find(opt.lambda > 0)
                switch i
                    case 2
                        fprintf('; WT: %6.3f',SP1f);
                    case 3
                        fprintf('; TV1: %6.3f',SP2f);
                    case 4
                        fprintf('; TV2: %6.3f',SP3f);
                    case 5
                        fprintf('; VS: %6.3f',VSf);
                    case 6
                        fprintf('; RI: %6.3f',RIf);
                    case 7
                        fprintf('; CD: %6.3f',CDf);
                end
            end
            fprintf('\n');
        end
        
        %   Update penalties to match data consistency target
        if opt.compression_levels > 0
            
            %   Compute average variance
            %DCf = DCf * FTnorm/numel(k);
            DCf = sqrt(sum(abs(knew(:)+t*kd(:) - k(:)).^2) / numel(k));
            
            %   Check if current compression level is close to target
            if (opt.compression_levels/DCf > 1.02 || opt.compression_levels/DCf<0.98)
                lmda = opt.lambda;
                opt.lambda = opt.lambda * (opt.compression_levels/DCf)^(0.4);
                grad_restart_flag = 1;
                
                fprintf('Relative variance: %3.2f\n',DCf/opt.compression_levels);
                for i = find(lmda>0)
                    fprintf('\tConstraint (%i): %7.6g -> %7.6g\n',i,lmda(i),opt.lambda(i));
                end
            end
            
        end
        

        %   Experimental
        if opt.lambda(6) && (nt >= opt.Npc) && (opt.RI.Iter/opt.update <= opt.Npc_Nup)
            if opt.verbose == 2
                fprintf('Perfomring principle component analysis\n');
            end
            %opt.IREF = svd_smooth(img,opt.Npc,0,opt.BF);
            opt.IREF = svd_smooth_patch(img,opt.Npc,0);
            grad_restart_flag = 1;
        end
        
        %   Reset conjugate gradient direction
        if grad_restart_flag
            grad_restart;
        end
        
        if opt.verbose == 2
            fprintf('-------\n');
        end
    end


%   Reset the conjugate gradient search direction when desired/needed
    function grad_restart
        delta = complex(zeros([np nv ns nt],opt.class),zeros([np nv ns nt],opt.class));
        if opt.verbose == 2
            fprintf('Conjugate Gradient Direction Reset\n');
        end
    end


%   Check the input variables for proper class
    function class_check
        
        %   Flag for printing to the screen
        updt = 0;
        
        %   Ensure single/float precision
        if strcmp(opt.class,'single')
            if ~isa(k,'single')
                k = single(k);
                updt = 1;
            end
            if ~isa(opt.vs_kern,'single')
                opt.vs_kern = single(opt.vs_kern);
            end
            if ~isa(opt.S,'single')
                opt.S = single(opt.S);
                updt = 1;
            end
            if ~isa(opt.kW,'single')
                opt.kW = single(opt.kW);
                updt = 1;
            end
            if ~isa(opt.IREF,'single')
                opt.IREF = single(opt.IREF);
                updt = 1;
            end
            
        %   Ensure double precision
        else
            if ~isa(k,'double')
                k = double(k);
                updt = 2;
            end
            if ~isa(opt.vs_kern,'double');
                opt.vs_kern = double(opt.vs_kern);
            end
            if ~isa(opt.S,'double');
                opt.S = double(opt.S);
                updt = 2;
            end
            if ~isa(opt.kW,'double')
                opt.kW = double(opt.kW);
                updt = 2;
            end
            if ~isa(opt.IREF,'double')
                opt.IREF = single(opt.IREF);
                updt = 2;
            end

        end
        
        %   Display info on the screen
        if opt.verbose
            if updt == 1
                fprintf('Inputs converted to single precision\n');
            elseif updt == 2
                fprintf('Inputs converted to double precision\n');
            end
        end
    end


%   Compute initial image
    function comp_img0
        if ~isfield(opt,'img0') || isempty(opt.img0)
            %   Condense all views if dynamic
            if nt > 1
                img = zeros(opt.size,opt.class);
                img(opt.U) = k;
                img = mean(img,4)./(mean(img~=0,4) + eps);
                img = opt.iFT(img,opt);
                img = repmat(img,[1 1 1 nt 1]);
                img = opt.iPH(opt.iPI(img,opt),opt);
                clear optt;
            else
                img = opt.iFTU(k,opt);
                img = opt.iPH(opt.iPI(img,opt),opt);
            end
        else
            img = opt.img0;
        end
        img(1) = img(1) + complex(0,eps); %  Ensures complex
        if strcmp(opt.class,'double') && ~isa(img,'double')
            img = double(img);
        elseif strcmp(opt.class,'single') && ~isa(img,'single')
            img = single(img);
        end
    end


%   Compute phase estimate
    function comp_Ph
        opt.phase_corr = opt.phase_corr(opt.phase_corr > 0);
        if any(opt.phase_corr)
            opt.Ph = opt.ePH(opt.fPH(img,opt),opt);
        end
    end

%   Check options structure against input data and disable certain options if warranted
    function param_check
        %   Check GRAPPA
        if  opt.lambda(1) > 0
            opt.lambda(1) = 0;
            if opt.verbose
                fprintf('GRAPPA constraint disabled\n');
            end
        end
        
        %   Check view-sharing
%         if nt == 1 && opt.lambda(5) && ~isfield(opt,'kvs')
%             opt.lambda(5) = 0;
%             if opt.verbose
%                 fprintf('View-sharing constraint disabled (only one time point)\n');
%             end
%         end
        
        %   Check that kspace weights exist
        if ~isfield(opt,'kW') || isempty(opt.kW)
            opt.kW = 1;
        end
        
        %   Check for reference image
        if opt.lambda(6) && isempty(opt.IREF)
            opt.IREF = complex(zeros([np nv ns nt],opt.class),zeros([np nv ns nt],opt.class));
        end
        
        %   Create variable to change conjugate gradient method
        opt.cg_mtd = 0;
    end


%   Estimate wavelet transform order
    function comp_worder
        if strcmp(opt.worder{1},'auto') && opt.lambda(2) > 0
            opt.worder = calc_wave_order(img(:,:,:,:,1),10);
            if opt.verbose
                fprintf('Estimated wavelet transform order: {');
                for i = 1:length(opt.worder)
                    fprintf('[%s]',num2str(opt.worder{i}));
                end
                fprintf('}\n');
            end
        end
    end


%   Compute l1-weighting factors
    function comp_l1weight
        if opt.l1weight > 0
            if opt.lambda(2)
                opt.WCwt  = estWT(opt.fSP1(img,opt.wname,opt.worder),opt);
            end
            if opt.lambda(3)
                opt.TV1wt = estWT(opt.fSP2(img,opt.FDorder1),opt);
            end
            if opt.lambda(4)
                opt.TV2wt = estWT(opt.fSP3(img,opt),opt);
            end
            if opt.lambda(5)
                opt.VSwt  = estWT(opt.fVS( img,opt),opt);
            end
            if opt.lambda(6)
                opt.RIwt  = estWT(opt.fSP4(img,opt),opt);
            end
            if opt.lambda(7)
                opt.CDwt  = estWT(opt.fCD(img,opt),opt);
            end
        else
            opt.WCwt  = [];
            opt.TV1wt = [];
            opt.TV2wt = [];
            opt.VSwt  = [];
            opt.RIwt  = [];
            opt.CDwt  = [];
        end
    end

end
