function [AIF, AIFROI, AIFparam, AIFtheta, AIFvp, vessels] = estimateAIF(img, opt)
%function [AIF, aifROI, param] = estimateAIF(img, opt)
%   Estimation of Vascular Input Functions   
%
%   Inputs:
%
%       img             (complex) image data [Nx Ny Nz Nt]
%       opt             options
%         .size         vector of dimension: [Nx Ny Nz Nt]
%         .FA           flip angle in degrees
%         .TR           TR in seconds
%         .TE           TE in seconds
%         .omega0       resonance frequency
%         .R1           R10 map [Nx Ny Nz]
%         .M0           M0 map [Nx Ny Nz]
%         .SIdirection  dimension of S/I direction
%         .ibolus       time index of bolus arrival
%         .time         time vector [Nt]
%         .voxel_size   voxel size in mm
%
%         .AIF.popAIF   population/site-specific AIF, set to [] if not used
%         .AIF.AIFROI   [x y z]-coordinates of ROI for AIF extraction,
%                       - each row is a new voxel
%                       - if set to [], will spawn ROI estimtion procedure
%         
%         .AIF.bloodT1  fixed blood T1, set to [] if not used
%         .AIF.bloodT2  fixed blood T2, set to [] if not used
%         .AIF.bloodr1  blood T1 relaxivity r1
%         .AIF.bloodr2  blood T2 relaxivity r2
%         .AIF.bloodchi blood susceptibility
%         .AIF.bloodHct Hematocrit value
%
%         .AIF.extractFrom          how to extract the AIF from image data:
%                                   'real'      image real part
%                                   'abs'       image magnitude
%                                   'phase'     image phase
%                                   'complex'   complex image data
%          .AIF.merge               how to compute the final AIF from the AIFs in
%                                   each voxel:
%                                   'mean'      just average all voxels
%          .AIF.paramterization     parameterization of the AIF
%                                   []                  no parameterized form
%                                   '2Gammas+Sigmoid'   
%
%   Define susceptibility options (for phase only)
%          .fSus = @(C, opt) (opt.chi/3)*C;     % conversion of concentration to susceptibility via dipole kernel
%          .iSus = @(C, opt) (3/(opt.chi))*C;   % this is an axis aligned dipole kernel, ie, only valid for vessels which are aligned with B0
%
%   Outputs:
%       AIF         probably what you are after when you call this function
%       AIFROI      linear indices of AIF ROI
%       AIFparam    (does not work yet) placeholder for the parameters if the
%                   AIF is forced into parameterized form
%       AIFtheta    map of angles between vessel and B0 (in degrees)
%       vessels     returns all vessels that have been detected
%
%   Example setup:
%       opt.fid                 = 1;                % log file handle
%       opt.size                = [256 256 160 50];
%       opt.FA                  = 15;
%       opt.TR                  = 0.005;
%       opt.TE                  = 0.002;
%       opt.SIdirection         = 1;
%       opt.ibolus              = 1;
%       opt.gamma               = 2*pi*42.58e6;     % gyromagnetic ratio 
%       opt.B0                  = 3;                % B0 field strength
%       opt.omega0              = opt.gamma*opt.B0;
%
%       opt.AIF.popAIF          = [];               % container for population AIF (empty set will try to estimate AIF from data)
%       opt.AIF.AIFROI          = [];               % keeping this empty will enforce ROI estimation
%       opt.AIF.extractFrom     = 'complex';        % 'abs' | 'real' | 'phase' | 'complex'
%       opt.AIF.merge           = 'mean';           % 'mean' | 'parameterized'
%       opt.AIF.paramterization = '2Gauss+Sigmoid'; % parameterization of the AIF
%       opt.AIF.integration     = 'sum';            % 'sum' | 'trapz' | 'adaptive'
%       opt.AIF.bloodT1         = [];               % [] (will use measured T1) | some fixed value that is used
%       opt.AIF.bloodT2         = 0.275;            % [] (will use measured T2 - currenlty not implemented) | some fixed value that is used
%           default values for blood T1 and T2star are given in Ref 3:
%               T1  = 1932ms;
%               T2* = 275ms;
%       opt.AIF.bloodr1         = 4.2;
%       opt.AIF.bloodr2         = 4.2;
%       opt.AIF.bloodchi        = 320e-9;           % molar susceptibility of the contrast agent: (roughly) 320 ppm/M = 320e-6 M^(-1) = 320e-9 L/mmol [2]
%       opt.AIF.bloodHct        = 0.45;             % average Hct value in US adult population according to Ref. [1]
%
%
%   References:
%       1. Bagher-Ebadian H, Jain R, Nejad-Davarani SP, et al. Model selection for DCE-T1 studies in glioblastoma. Magn. Reson. Med. 2012;68:241?251 doi: 10.1002/mrm.23211.
%       2. J Korporaal, et al. "Phase-based arterial input function measurements in the femoral arteries for quantification of dynamic contrast-enhanced (DCE) MRI and comparison with DCE-CT," MRM, 2011, 1267-1274.
%       3. Klawer EME, van Houdt PJ, Simonis FFJ, et al. Improved repeatability of dynamic contrast-enhanced MRI using the complex MRI signal to derive arterial input functions: a test-retest study in prostate cancer patients. Magn. Reson. Med. 2019:1?12 doi: 10.1002/mrm.27646.
%       4. Fluckiger JU, Schabel MC, DiBella EVR. Model-based blind estimation of kinetic parameters in Dynamic Contrast Enhanced (DCE)-MRI. Magn. Reson. Med. 2009;62:1477?1486.
%       5. Schabel MC, Fluckiger JU, DiBella EVR, et al. A model-constrained Monte Carlo method for blind arterial input function estimation in dynamic contrast-enhanced MRI: I. simulations. Phys. Med. Biol. 2010;55:4783?4806 doi: 10.1088/0031-9155/55/16/011.
%       6. Schabel MC, DiBella EVR, Jensen RL, Salzman KL. A model-constrained Monte Carlo method for blind arterial input function estimation in dynamic contrast-enhanced MRI: II. In vivo results. Phys. Med. Biol. 2010;55:4807?4823 doi: 10.1088/0031-9155/55/16/012.
%       7. Ji?ík R, Taxt T, Mací?ek O, et al. Blind deconvolution estimation of an arterial input function for small animal DCE-MRI. Magn. Reson. Imaging 2019;62:46?56 doi: 10.1016/j.mri.2019.05.024.
%
% Yannick Bliesener 2019
%

AIFparam = [];
AIFtheta = [];
vessels  = [];
AIFIdx   = [];

% Determine AIF ROI
if ~isempty(opt.AIF.popAIF)
    AIF = opt.AIF.popAIF;
    AIFROI = [];
    return
elseif isempty(opt.AIF.AIFROI)
    [AIFROI, vessels, AIFIdx] = estimate_AIF_ROI(img, opt);
else
    AIFROI = opt.AIF.AIFROI;
end

if isempty(AIFROI)
    err_msg = 'Could not find AIF ROI!';
    fprintf(opt.fid, '\t\t\t ERROR: %s \n', err_msg);
    
    AIF = [];
    AIFROI = [];
    return
end

% get potential AIFs
img = reshape(img, prod(opt.size(1:3)), opt.size(4));
AIFsignal = img(AIFROI(:), :);

% copy parameters
optA.size    = [size(AIFsignal,1), 1, 1, size(AIFsignal,2)];
optA.FA      = opt.FA;
optA.TR      = opt.TR;
optA.TE      = opt.TE;
optA.omega0  = opt.omega0;
optA.TK.r1   = opt.AIF.bloodr1;
optA.TK.r2   = opt.AIF.bloodr2;
optA.TK.chi  = opt.AIF.bloodchi;
optA.bloodvp = 1-opt.AIF.bloodHct;

% B1
if isfield(opt,'B1') && ~isempty(opt.B1)
    optA.B1 = opt.B1(AIFROI(:));
end

% get blood T1
if isempty(opt.AIF.bloodT1)
    optA.R1 = opt.R1(AIFROI(:));
    optA.M0 = opt.M0(AIFROI(:));
else
    optA.R1 = 1./opt.AIF.bloodT1 .* ones(numel(AIFROI),1);
    
    % estimate M0
    baseline = mean(AIFsignal(:,1:opt.ibolus), 2);
    optA.M0  = abs(baseline) ./ SPGR(1, optA.R1, optA);
end

% get blood T2
if isempty(opt.AIF.bloodT2)
    optA.R2 = opt.R2(AIFROI(:));
else
    optA.R2 = 1./opt.AIF.bloodT2 .* ones(numel(AIFROI),1);
end

% Theta (angle of vessel with B0 in degrees)
if isfield(opt,'theta') && ~isempty(opt.theta)
    optA.theta = opt.theta(AIFROI(:));
else
    optA.theta = [];
end

% tissue TK if given
if isfield(opt.TK, 'vp') && ~isempty(opt.TK.vp)
    optA.TK.vp = opt.TK.vp(AIFROI(:));
    optA.TK.kt = opt.TK.kt(AIFROI(:));
    optA.TK.kep = opt.TK.kep(AIFROI(:));
end

% extract AIFs
switch opt.AIF.extractFrom
        
    case {'real', 'abs'}
        if strcmp(opt.AIF.extractFrom, 'real')
            AIFsignal = real(AIFsignal);
        else
            AIFsignal = abs(AIFsignal);
        end
        
        % estimate the AIF from all the given tissue parameters
        switch opt.AIF.AIFROItype
            case 'tissue'
                % this approach uses current estimates of TK parameters to
                % estimate the AIF from all enhancing tissue, see Ref. 4-7.
                % there are probably many more references and flavors in
                % the regard, essentially the idea is that blood has too
                % much iron, hence susceptibilty, and too much partial
                % volume averaging, so estimating the AIF from other
                % enhancing tissues might be a good idea.
                
                conc = iSPGR(AIFsignal, optA);
                indError = sum(isnan(conc) | isinf(conc), 2) > 0;
                conc = iFXL(conc, optA);
                
                % erase model errors from the AIF estimation
                conc(indError, :) = [];
                
                optA.tol     = 1e-4;
                optA.maxIter = 200;
                optA.solver  = 'pcg';
                optA.init    = opt.AIF.initAIF;
                optA.size    = [size(conc,1), 1, 1, size(conc,2)];
                optA.time    = opt.time;
                
                optA.TK.vp(indError)  = [];
                optA.TK.kt(indError)  = [];
                optA.TK.kep(indError) = [];

                AIF = estimateAIFfromTissueConcentration(conc, optA);
            case 'vessel'
                optA.lambda  = opt.AIF.lambda;
                    % the lambda parameter controls the allowed deviation of
                    % vessel vps from the nominal blood vp
                
                [AIF, AIFvp] = estimateAIFfromMagnitude(AIFsignal, optA);
%                 AIF = mean(conc, 1);
%                 AIF = AIF ./ opt.AIF.bloodvp;
            otherwise
                err_msg = sprintf('Unkown reference tissue for %s AIF extraction!', opt.AIF.extractFrom);
                fprintf(opt.fid, '\t\t\t ERROR: %s \n', err_msg);
                error(err_msg)
        end
        
    case 'phase'
        % Estimate the angle between vessel and B0
        if isempty(optA.theta)
            if isempty(AIFIdx)
                optA.theta = estimateVesselTheta(AIFROI, [], opt);
            else
                theta = estimateVesselTheta(vessels{AIFIdx}, AIFROI, opt);
                tmp_theta = nan(opt.size(1:3));
                tmp_theta(AIFROI(:)) = theta;
                optA.theta = tmp_theta(AIFROI(:));
            end
        end
        
        AIFtheta = nan(opt.size(1:3));
        AIFtheta(AIFROI(:)) = optA.theta;
        
        [~, AIFs] = iSPGR(AIFsignal, optA);
        
        % dipole kernel
        F = reshape(1/3 - 0.5*sin(optA.theta*pi/180).^2, [size(AIFs,1) 1]);
        
        AIF = F \ AIFs;
        AIF = AIF ./ optA.TK.chi;
        
        % estimate the AIF from all the given tissue parameters
        switch opt.AIF.AIFROItype
            case 'tissue'
                error('Estimating vascular input funtions from the phase of all enhancing tissue requires implementation! :/')
            case 'vessel'
                AIF = AIF ./ optA.bloodvp;
                AIFvp = optA.bloodvp;
            otherwise
                err_msg = sprintf('Unkown reference tissue for %s AIF extraction!', opt.AIF.extractFrom);
                fprintf(opt.fid, '\t\t\t ERROR: %s \n', err_msg);
                error(err_msg)
        end
        
    case 'complex'
        
        % Estimate the angle between vessel and B0
        if isempty(optA.theta)
            if isempty(AIFIdx)
                optA.theta = estimateVesselTheta(AIFROI, [], opt);
            else
                theta = estimateVesselTheta(vessels{AIFIdx}, AIFROI, opt);
                tmp_theta = nan(opt.size(1:3));
                tmp_theta(AIFROI(:)) = theta;
                optA.theta = tmp_theta(AIFROI(:));
            end
        end
        
        AIFtheta = nan(opt.size(1:3));
        AIFtheta(AIFROI(:)) = optA.theta;
        
        % estimate the AIF from all the given tissue parameters
        switch opt.AIF.AIFROItype
            case 'tissue'
                error('Estimating vascular input funtions from all enhancing tissue in complex plane requires implementation! :/')
            case 'vessel'
                optA.lambda  = opt.AIF.lambda;
                    % the lambda parameter controls the allowed deviation of
                    % vessel vps from the nominal blood vp
                
                [AIF, AIFvp] = estimateAIFinComplexPlane(AIFsignal, optA);
            otherwise
                err_msg = sprintf('Unkown reference tissue for %s AIF extraction!', opt.AIF.extractFrom);
                fprintf(opt.fid, '\t\t\t ERROR: %s \n', err_msg);
                error(err_msg)
        end
        
    otherwise
        err_msg = 'Unkown data type to extract AIF!';
        fprintf(opt.fid, '\t\t\t ERROR: %s \n', err_msg);
        error(err_msg)
end

% fit to AIF model
if ~isempty(opt.AIF.paramterization)
    [ AIFparam ] = parameterizeAIF( AIF, opt );
end

end