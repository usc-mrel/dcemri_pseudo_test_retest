function [TK,img,outputInfo] = MOCCOTK(k,opt)
%function [TK,img,outputInfo] = MOCCOTK(k,opt)
%   plain vanilla MOCCO TK estimation
%
%   Inputs:
%       k       k-space data
%       opt     reconstruction options, see TKESTIMATION_optset.m
%
%   Outputs:
%       TK          container of tracer kinetic parameter maps
%       img         estimated image time series
%       outputInfo  auxiliary information
%
%   References:
%       1. Guo Y, Lingala SG, Bliesener Y, Lebel RM, Zhu Y, Nayak KS. Joint arterial input function and tracer kinetic parameter estimation from undersampled dynamic contrast-enhanced MRI using a model consistency constraint. Magn. Reson. Med. 2018;79:2804?2815 doi: 10.1002/mrm.26904.
%
%   Authors:
%       Yannick Bliesener, bliesene@usc.edu, MREL at USC, 2019
%

%% Setup
Nx = opt.size(1);
Ny = opt.size(2);
Nz = opt.size(3);
Nt = opt.size(4);
Nc = opt.size(5);
Nk = numel(k);

outputInfo.DCcost    = -ones(opt.MaxIter,1);
outputInfo.MOCCOcost = -ones(opt.MaxIter,1);
outputInfo.cost      = -ones(opt.MaxIter,1);
outputInfo.AIFs      = -ones(opt.MaxIter,Nt);
outputInfo.AIFparams = cell(opt.MaxIter,1);
outputInfo.AIFROI    = cell(opt.MaxIter,1);
outputInfo.AIFtheta  = cell(opt.MaxIter,1);
outputInfo.tforms    = cell(opt.MaxIter,1);
outputInfo.Vessels   = cell(opt.MaxIter,1);

% look for GPUs if available
if (gpuDeviceCount > 0)
    opt.use_gpu = 1;
    %     parpool(gpuDeviceCount);
    % matlab by default puts one GPU per worker
    
    opt.chunkSize = 20;
    % larger is faster but takes more memory
    % especially the SENSE recon is craving for memory for some reason
    opt.num_of_chunks = ceil(Nx/opt.chunkSize);
    
    gpu_id = opt.gpuid;
    gpuDevice(gpu_id);
else
    opt.use_gpu = 0;
end

% if exist('debug.mat', 'file')
%     fprintf(opt.fid,'\tLoading debugging information from file....\n');
%     load('debug.mat', 'img')
% else
    fprintf(opt.fid,'\tInitial MOCCO POCS iterations...');tic
    if opt.use_gpu
        % all this drama just to make the data fit on GPUs
        % in case you have a GPU with mem > 16GB onemight be able to relax this
        
        % if you have several GPUs you could switch the for loop to a
        % parfor loop and this will parallize over GPUs
        
        utmp = false(opt.size);
        utmp(opt.U) = 1;
        
        img = cell(opt.num_of_chunks,1);
        change = cell(opt.num_of_chunks,1);
        for ch=1:opt.num_of_chunks
            chunkInd = 1+(ch-1)*opt.chunkSize:min(ch*opt.chunkSize,Nx);
            
            gk             = gpuArray(k(chunkInd,:,:,:,:));
            
            gOpt = struct();
            gOpt.SENSE.kspaceDataType = opt.SENSE.kspaceDataType;
            gOpt.fPI  = opt.fPI;
            gOpt.iPI  = opt.iPI;
            gOpt.fFT  = opt.fFT;
            gOpt.iFT  = opt.iFT;
            gOpt.iFTU = opt.iFTU;
            gOpt.FTshift = opt.FTshift;
            gOpt.FTdim = opt.FTdim;
            gOpt.use_gpu = opt.use_gpu;
            
            gOpt.brainmask = gpuArray(opt.brainmask(chunkInd,:,:));
            gOpt.U         = gpuArray(find(utmp(chunkInd,:,:,:,:)));
            gOpt.S         = gpuArray(opt.S(chunkInd,:,:,:,:));
            gOpt.size      = [length(chunkInd), opt.size(2:end)];
            
            [tmp_img, tmp_change] = initPOCS(gk, gOpt);
            
            img{ch} = gather(tmp_img);
            change{ch} = gather(tmp_change);
            
        end
        img = cell2mat(img);
        
        clear gOpt
        gpuDevice(gpu_id);
        wait(gpuDevice);
        % wait until GPU is cleared before sending new data
        % this should get rid of the memory spikes between iterations
        % and can result in bigger chunk sizes being used.
        
        clear gk gOpt
        
%         save('debug.mat', 'img', 'opt', 'change', '-v7.3')
    else
        img = initPOCS(k, opt);
    end
    fprintf(opt.fid,'done (%g sec)\n',toc);
% end
% baseline = SPGR(opt.M0, opt.R1, opt);

% copy options for SENSE
optSENSE = opt;
optSENSE.solver         = opt.SENSE.solver;
optSENSE.MaxIter        = opt.SENSE.MaxIter;
optSENSE.tol            = opt.SENSE.tol;
optSENSE.kspaceDataType = opt.SENSE.kspaceDataType;
% lambda is set during the iteration

% setup motion correction
tform  = [];
tformI = [];

% setup containers
TK        = [];
modelMask = [];

%%  MOCCO iterations
Dobj = ones(1,opt.Dobjn);

fprintf(opt.fid,'\tCore MOCCO iterations:\n\t');tic
for i = 1:opt.MaxIter
    
    %% Set lambda
    if length(opt.MOCCO.lambda) > 1
        optSENSE.lambda = opt.MOCCO.lambda(i);
    else
        optSENSE.lambda = opt.MOCCO.lambda;
    end
    
    %%   Re-estimate phase correction term
    % this should make the AIF estiamtion more stable and also help with
    % CTC conversion
    
    opt.Ph = repmat(exp(sqrt(-1)*angle(img(:,:,:,1))), [1 1 1 opt.size(4)]);
    img = opt.iPH(img,opt);
    
    %% MOCCO: AIF estimation
    aifOpt = opt;
    aifOpt.S = [];
    aifOpt.U = [];
    if (i > 1) && strcmp(aifOpt.AIF.AIFROItype, 'wholebrain')
        aifOpt.TK.vp = TK.vp;
        aifOpt.TK.kt = TK.kt;
        aifOpt.TK.kep = TK.kep;
        
        [~,ibolus] = min(abs(aifOpt.AIF.shifts));
        aifOpt.AIF.AIFROI = aifOpt.brainmask & ...
            (TK.shift == ibolus) & ...                      % the code cannot handle bolus delay
            ((TK.vp > 0) | (TK.kt > 0)) & ...               % enhancing tissue
            (TK.vp < aifOpt.TK.parameterBounds.vp(2)) & ...    % exlcude model fitting errors
            (TK.kt < aifOpt.TK.parameterBounds.kt(2));         % exlcude model fitting errors
            
        aifOpt.AIF.initAIF = AIF;

    end
    
%     if (i == 1) || (i == 2)
%         aifOpt.AIF.extractFrom = 'abs';
%     end
    
    [AIF, AIFROI, AIFparam, AIFtheta, AIFvp, vessels] = estimateAIF(img, aifOpt);
    
    outputInfo.AIFROI{i} = AIFROI;
    if ~isempty(AIF)
        outputInfo.AIFs(i,:) = AIF;
        outputInfo.AIFparams{i} = AIFparam;
    end
    if ~isempty(AIFtheta)
        outputInfo.AIFtheta{i} = AIFtheta(AIFROI(:));
    end
    
    % replace AIF if parameterization is available
    if ~isempty(opt.AIF.paramterization)
        switch opt.AIF.paramterization
            case '2Gammas+Sigmoid'
                fittedAIF = AIFgammavariates( opt.time, AIFparam );
            otherwise
                err_msg = 'MOCCOTK: Unkown AIF model to fit to!';
                fprintf(opt.fid, '\t ERROR: %s \n', err_msg);
                error(err_msg)
        end
        
        % Check if fit worked out well
        [FMI, FRI] = balvaysCriterion(AIF,fittedAIF);
        if (FMI > 0.99)
            % this means that the fit should be able to explain more than 99% of the information in the measurement
            AIF = fittedAIF;
            
            if i >= 2
                % use a good fit as initial value for the next round
                opt.AIF.initialParameters = AIFparam;
            end
        else
            fprintf(opt.fid, '\tWarning: %s \n', 'Fit to AIF has been rejected!');
        end
        
        tmask = ( AIF > eps );
        Nt = sum( tmask );
        
        if Nt < 2
            err_msg = 'MOCCOTK: AIF fitting failed!';
            fprintf(opt.fid, '\t ERROR: %s \n', err_msg);
            break
        end
    end
    
    if isempty(AIF) || any(isnan(AIF)) || any(isinf(AIF))
        fprintf(opt.fid,'\tFailed to estimate AIF\n');
        break
    end
    
    figure(1), plot( opt.time, AIF );
    title('Current AIF estimate')
    xlabel('Time')
    drawnow
    
    if i == 3
        % keep this for the rest of the iterations
        opt.AIF.AIFROI = AIFROI;
        opt.theta = AIFtheta;
    end
    
    %% Retrieve concentration time curves
    img = abs(img);
    
    % inter frame motion correction
    if opt.MC.interframe && (i==20)
        %         [~,tform] = registervolumes(img ,opt.FOV, []);
        %         T = -ones(4,4,size(img,4));
        %         for t=1:size(img,4)
        %             T(:,:,t) = tform(t).T;
        %         end
        %         T = reshape(T, 16, size(img,4));
        %         for t=1:16
        %             T(t,:) = medfilt1(T(t,:), [], 2);
        %         end
        %         T = reshape(T, [4, 4, size(img,4)]);
        %         for t=1:size(img,4)
        %             T(:,4,t) = [0 0 0 1];
        %             tform(t).T = T(:,:,t);
        %         end
        
    end
    %     if mod(i,2) == 0
    % re-estimate motion
    %         tform = [];
    %     end
    
    %     if opt.MC.interframe && (i>=10)
    %         [img,tform,tformI,Mot] = registervolumes(img ,opt.FOV, tform, [], 'monomodal');
    %     end
    outputInfo.tforms{i} = tform;
    
    % Signal to R1 map
    R1 = iSPGR(img, opt);
    
    % R1 to CTC
    conc = iFXL(R1, opt);
    conc( isinf(conc) | isnan(conc) ) = 0;
    
    %% MOCCO: TK estimation
    
    if (i == 1) || (i == 2)
        % will do Patlak/LLSQ init for ETK
        initial_param = [];
    else
        initial_param.kt    = TK.kt;
        initial_param.vp    = TK.vp;
        initial_param.ve    = TK.ve;
        initial_param.kep   = TK.kep;
        initial_param.fp    = TK.fp;
        initial_param.shift = TK.shift;
    end
    
    [TK, predSignal, OUT] = estimateTK(conc, AIF, initial_param, opt);
    modelMask = OUT.modelMask;
    
    if gpuDeviceCount > 0
        % clear all stuff from GpuFit
        gpuDevice(gpu_id);
        wait(gpuDevice);
    end
    
    %% MOCCO: compute predicted signal
    
    if opt.MOCCO.enforceVesselParam
        
        AIFvp(isnan(AIFvp)) = 1 - opt.AIF.bloodHct;
        
        TK.vp(AIFROI(:)) = AIFvp(:) * (opt.AIF.bloodr1 / opt.TK.r1);
        TK.kt(AIFROI(:)) = 0;
        TK.ve(AIFROI(:)) = 1 - TK.vp(AIFROI(:));
        [~,ibolus] = min(abs(opt.AIF.shifts));
        TK.shift(AIFROI(:)) = ibolus;
        
        try
            if isempty(vessels)
                [~, vessels] = estimate_AIF_ROI(img, opt);
            end
            outputInfo.Vessels{i} = vessels;
        catch e
            outputInfo.Vessels{i} = [];
        end
        
        predSignal = tk2conc(TK, AIF, opt);
    end
    
    predR1 = fFXL(predSignal, opt);
    if ~isempty(opt.fSus)
        deltaB = opt.fSus(predSignal, opt);
    else
        deltaB = [];
    end
    predSignal = fSPGR(predR1, [], deltaB, opt);
    % R2 is unused
    
    % Re-apply motion correction
    %     if opt.MC.interframe && (i>=10)
    %         predSignal = registervolumes(predSignal,opt.FOV,tformI);
    %         img = registervolumes(img,opt.FOV,tformI);
    %     end
    
    %   Re-apply phase correction term
    predSignal = opt.fPH(predSignal,opt);
    img = opt.fPH(img,opt); % used as initial guess
    
    %% MOCCO: SENSE part
    
    if opt.use_gpu
        
        curDCcost = 0;
        
        utmp = false(opt.size);
        utmp(opt.U) = 1;
        
        imgCell = cell(opt.num_of_chunks,1);
        for ch=1:opt.num_of_chunks
            chunkInd = 1+(ch-1)*opt.chunkSize:min(ch*opt.chunkSize,opt.size(1));
            
            gk             = gpuArray(k(chunkInd,:,:,:,:));
            gPS            = gpuArray(predSignal(chunkInd,:,:,:,:));
            gInit          = gpuArray(img(chunkInd,:,:,:,:));
            
            gOpt = struct();
            gOpt.SENSE.kspaceDataType = optSENSE.kspaceDataType;
            gOpt.fPI  = opt.fPI;
            gOpt.iPI  = opt.iPI;
            gOpt.fFT  = opt.fFT;
            gOpt.iFT  = opt.iFT;
            gOpt.iFTU = opt.iFTU;
            gOpt.fFTU = opt.fFTU;
            gOpt.FTshift = opt.FTshift;
            gOpt.FTdim = opt.FTdim;
            
            gOpt.solver  = optSENSE.solver;
            gOpt.MaxIter = optSENSE.MaxIter;
            gOpt.tol     = optSENSE.tol;
            gOpt.lambda  = optSENSE.lambda;
            
            gOpt.U         = gpuArray(find(utmp(chunkInd,:,:,:,:)));
            gOpt.S         = gpuArray(opt.S(chunkInd,:,:,:,:));
            gOpt.size      = [length(chunkInd), opt.size(2:end)];
            
            [gImg,OUT] = MOCCO_recon(gk,gPS,gInit,gOpt);
            
            imgCell{ch} = gather(gImg);
            
            %   Compute data consistency on GPU
            gImg = fSENSE(gImg,gOpt);
            gImg = fFFTU(gImg,gOpt);
            tmpDC = sum(abs(gImg(:) - gk(:)).^2);
            
            curDCcost = curDCcost + gather(tmpDC);
            
        end
        gpuDevice(gpu_id);
        wait(gpuDevice);
        % wait until GPU is cleared before sending new data
        
        img = cell2mat(imgCell);
        clear imgCell
        
    else
        [img,OUT] = MOCCO_recon(k,predSignal,img,optSENSE);
        
        %   Compute data consistency
        gInit = fSENSE(img,opt);
        gInit = fFFTU(gInit,opt);
        
        curDCcost = sum(abs(gInit(:) - k(:)).^2);
        
    end
    
    
    %%  Plotting
    figure(2);
    subplot(1,2,1);
    vpt = squeeze(TK.vp(ceil(size(TK.vp,1)/2),:,:));
    imagesc(vpt,[0 0.6]);
    colormap jet;axis square off;colorbar;
    title('vp map')
    subplot(1,2,2);
    ktt = squeeze(TK.kt(ceil(size(TK.kt,1)/2),:,:)) * 60;
    imagesc(ktt,[0 0.6]);
    colormap jet;axis square off;colorbar;
    title('kt map')
    
    imgt = abs(squeeze(img(ceil(size(img,1)/2),:,:,:)));
    imgt = reshape(imgt,[size(imgt,1) size(imgt,2) 1 size(imgt,3)]);
    imgt = 0.95*imgt./max(imgt(:));
    imgt = repmat(imgt,[1 1 3 1]);
    figure(3);
    montage(imgt);
    drawnow;
    title('Temporal evolution')
    
    if ~isempty(AIFROI)
        imgt = abs(squeeze(img(:,:,:,opt.ibolus+2)));
        aifroit = false(opt.size(1:3));
        aifroit(AIFROI(:)) = 1;
        figure(4);
        montage_with_ROI(permute(imgt, [2 3 1]), permute(aifroit, [2 3 1]), []);
        drawnow;
        title('AIF ROI')
    end
    
    %%  compute costs
    
    outputInfo.MOCCOcost(i) = norm( predSignal(:) - img(:), 2).^2;
    clear predSignal
    
    %   compute cost
    outputInfo.DCcost(i) = curDCcost;
    outputInfo.cost(i)  = outputInfo.DCcost(i) + optSENSE.lambda*outputInfo.MOCCOcost(i);
    
    % 	if i > 1
    %         fprintf('after SENSE cost: %e\n', outputInfo.cost(i))
    %     end
    
    
    %   Track changes in the objective function
    if i > 1
        Dobj = circshift(Dobj,[0 -1]);
        Dobj(end) = 100*(outputInfo.cost(i-1)-outputInfo.cost(i))./outputInfo.cost(i);
        
        if (i > opt.MinIter) && median(Dobj) < opt.Dobj
            fprintf(opt.fid,'\t\tInsufficient progress in iteration %i!\n',i);
            break
        end
    end
    
    % Plotting
    figure(5);
    subplot(3,1,1);
    semilogy(1:i,outputInfo.DCcost(1:i)./outputInfo.DCcost(1),'x-');
    grid on;
    drawnow;
    title('Data consistency')
    subplot(3,1,2);
    semilogy(1:i,outputInfo.MOCCOcost(1:i)./outputInfo.MOCCOcost(1),'x-');
    grid on;
    drawnow;
    title('Model consistency')
    subplot(3,1,3);
    semilogy(1:i,outputInfo.cost(1:i)./outputInfo.cost(1),'x-');
    grid on;
    drawnow;
    title('Total cost')
    
    fprintf(opt.fid,' %i', i);
    
    fprintf('\n')
    
end
fprintf(opt.fid,'\n\tdone (%g sec)\n',toc);

outputInfo.iterations = i;

%truncate output
outputInfo.cost         = outputInfo.cost(1:i);
outputInfo.MOCCOcost    = outputInfo.MOCCOcost(1:i);
outputInfo.DCcost       = outputInfo.DCcost(1:i);
outputInfo.AIFs         = outputInfo.AIFs(1:i,:);
outputInfo.AIFparams    = outputInfo.AIFparams(1:i);
outputInfo.AIFROI       = outputInfo.AIFROI(1:i);
outputInfo.AIFtheta     = outputInfo.AIFtheta(1:i);
outputInfo.modelMask    = modelMask;
outputInfo.tforms       = outputInfo.tforms(1:i);
outputInfo.tform        = tform;
outputInfo.tformI       = tformI;
outputInfo.Vessels      = outputInfo.Vessels(1:i);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Misc nested functions %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [img, change] = initPOCS(k, opt)

% initial guess
switch opt.SENSE.kspaceDataType
    case 'full'
        img = opt.iFT(k,opt);
    case 'list'
        img = opt.iFTU(k,opt);
end
img = opt.iPI(img,opt);

num_iter = 50;
change = -ones(num_iter,1);
if opt.use_gpu
   change = gpuArray(change); 
end
% initial POCS to get some starting AIF
for ii=1:num_iter
    
    prev_img = img;
    
    % apply brain mask
    img = reshape(img, prod(opt.size(1:3)), opt.size(4));
    img( ~opt.brainmask(:), : ) = 0;
    img = reshape(img, opt.size(1:4));
    
    % improve initial guess
    img = reshape(img, [], opt.size(4));
    img = truncateSingValMatrix(img, 8);
    img = reshape(img, opt.size(1:4));
    
    %   Convert to k-space
    img = opt.fPI(img,opt);
    img = opt.fFT(img,opt);
    
    % replace samples
    img(opt.U) = k;
    
    %   Convert to coil-combined image
    img = opt.iFT(img,opt);
    img = opt.iPI(img,opt);
    
    change(ii) = norm(prev_img(:) - img(:));
    
end
end

