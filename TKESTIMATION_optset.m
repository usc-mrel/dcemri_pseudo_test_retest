function opt = TKESTIMATION_optset
%   These are the default parameter settings for TK parameter estimation
%   with the STARDCE package
%
%
%
%
%   References:
%        1. Bagher-Ebadian H, Jain R, Nejad-Davarani SP, et al. Model selection for DCE-T1 studies in glioblastoma. Magn. Reson. Med. 2012;68:241?251 doi: 10.1002/mrm.23211.
%        2. J Korporaal, et al. "Phase-based arterial input function measurements in the femoral arteries for quantification of dynamic contrast-enhanced (DCE) MRI and comparison with DCE-CT," MRM, 2011, 1267-1274.
%        3. Klawer EME, van Houdt PJ, Simonis FFJ, et al. Improved repeatability of dynamic contrast-enhanced MRI using the complex MRI signal to derive arterial input functions: a test-retest study in prostate cancer patients. Magn. Reson. Med. 2019:1?12 doi: 10.1002/mrm.27646.
%        4. Pintaske J, Martirosian P, Graf H, et al. Relaxivity of Gadopentetate Dimeglumine (Magnevist), Gadobutrol (Gadovist), and Gadobenate Dimeglumine (MultiHance) in Human Blood Plasma at 0.2, 1.5, and 3 Tesla. Invest. Radiol. 2006;41:213?221 doi: 10.1097/01.rli.0000197668.44926.f7.
%        5. Simonis FFJ, Sbrizzi A, Beld E, Lagendijk JJW, van den Berg CAT. Improving the arterial input function in dynamic contrast enhanced MRI by fitting the signal in the complex plane. Magn. Reson. Med. 2016;76:1236?1245 doi: 10.1002/mrm.26023.
%        6. Sourbron SP, Buckley DL. Classic models for dynamic contrast-enhanced MRI. NMR Biomed. 2013;26(8):1004-1027.
%        7. Tofts P. T1-weighted DCE Imaging Concepts: Modelling, Acquisition and Analysis. Signal. 2010;500(450):400.
%
%   Yannick Bliesener 2019


%%   Define general options
opt.MaxIter = 100;              %   Maximum iterations
opt.MinIter = 30;               %   Minimum iterations
opt.class   = 'single';         %   Data class for computation (double or single)
opt.gpu_id  = [];               %   If multiple GPUs are available this can be used to pin a specific GPU by id
opt.Dobj    = 1e-3;             %   Stopping threshold (% chanage in objective)
opt.Dobjn   = 6;                %   Median window for changes below threshold
opt.plot    = 0;                %   Flag for plotting
opt.plotf   = @(img,opt)recon_plot(img,opt);  %   Plotting function
opt.verbose = 2;                %   Flag for text output (0: none; 1: some; 2: lots)
opt.update  = 5;                %   Periodic update interval
opt.size    = [];               %   Image size (placeholder)
opt.U       = [];               %   Sampled indices in k-space

%% Global constants

opt.gamma        = 2*pi*42.58e6;        % gyromagnetic ratio 
opt.B0           = 3;                   % B0 field strength
opt.omega0       = opt.gamma*opt.B0;
opt.SIdirection  = 1;                   % dimension of S/I direction


%% Operators

%   Define phase correction options
opt.phase_corr = [];                    %   Phase correction dimensions
opt.Ph_cf = [0.25 0.05 0.05];           %   Cutoff frequencies for phase estimation
opt.ePH = @(img,opt)estPH(img,opt);     %   Phase estimation function
opt.fPH = @(img,opt)fPH(img,opt);       %   Forward phase correction
opt.iPH = @(img,opt)iPH(img,opt);       %   Inverse phase correction

%   Define FT options
opt.fFT = @(f,opt)fFastFT(f,opt);       %   Define the forward Fourier operator
opt.iFT = @(f,opt)iFastFT(f,opt);       %   Define the inverse Fourier operator
opt.fFTU= @(f,opt)fFFTU(f,opt);         %   Define the forward Fourier operator (undersampled)
opt.iFTU= @(f,opt)iFFTU(f,opt);         %   Define the inverse Fourier operator (undersampled)
opt.FTshift = 0;                        %   FT shift
opt.FTdim = [2 3];                      %   FT dimension

%   Define GRAPPA/SPIRiT options
opt.fPI = @(img,opt)fSENSE(img,opt);    %   Forward SENSE transform
opt.iPI = @(img,opt)iSENSE(img,opt);    %   Inverse SENSE transform

%   Define susceptibility options
opt.fSus = @(C,opt) fSus(C, opt);       % conversion of concentration to susceptibility via dipole kernel
opt.iSus = [];                          % this is an axis aligned dipole kernel, ie, only valid for vessels which are aligned with B0

%% Parameters for various submodules

%   Define TK estimation settings
opt.TK.models               = {'patlak', 'etk'};                % 'null' | 'vpmodel' | 'patlak' | 'tofts'| 'etk' | 'tcxm'
                                                                % NOTE:
                                                                %   - setting multiple models will enforce model selection
                                                                %   - if multiple models are used it is very important
                                                                %     that the most complex model comes last in the list!!
opt.TK.selectionCriterion   = 'AIC';                            % selection criterion for model selection: 'AIC' | 'BIC' | 'Balvays' | 'SSE'
opt.TK.r1                   = 4.5;                              % relaxivity [7]
opt.TK.chi                  = 320e-9;                           % susceptibility of contrast agent [2]
opt.TK.parameterBounds.vp   = [0 1];                            % [min max]
opt.TK.parameterBounds.ve   = [0.02 1];
opt.TK.parameterBounds.kt   = [0 inf];                         % in 1/s
opt.TK.parameterBounds.kep  = [0 inf];                          % in 1/s
opt.TK.parameterBounds.fp   = [0.001 100/60];                   % in ml/s/100ml (this is different from the usual: ml/min/100ml)

opt.TK.tofts.initialFit     = 'LLSQ';                           % determines the initialization for non-linear fitting of tofts:
                                                                % 'patlak' | 'LLSQ'
                                                                % 'LLSQ' is linear fitting as described in Murase MRM 2004 and Flouri MRM 2016
opt.TK.tofts.solver         = 'gpufit';                         % fitting routine for Tofts model
                                                                % 'gpufit' | 'LLSQ' | 'ncg'

opt.TK.etk.initialFit       = 'LLSQ';                           % determines the initialization for non-linear fitting of etk/tcxm:
                                                                % 'patlak' | 'LLSQ'
                                                                % 'LLSQ' is linear fitting as described in Murase MRM 2004 and Flouri MRM 2016
opt.TK.etk.solver           = 'gpufit';                         % fitting routine for ETK model
                                                                % 'gpufit' | 'LLSQ' | 'ncg' |'gaussnewton'
opt.TK.etk.init.kep         = [];                               % [] | array of values for kep in 1/s
                                                                % if not empty overwrites the value from opt.TK.etk.initialFit for kep and
                                                                % refits to retain best fit
                                                                % Example:
                                                                % opt.TK.etk.init.kep = linspace(0, 1/60, 10);
                                                    
opt.TK.tcxm.initialFit      = 'LLSQ';                           % determines the initialization for non-linear fitting of etk/tcxm:
                                                                % 'LLSQ'
                                                                % 'LLSQ' is linear fitting as described in Murase MRM 2004 and Flouri MRM 2016
opt.TK.tcxm.solver          = 'gpufit';                         % fitting routine for 2CXM model
                                                                % 'gpufit' | 'LLSQ'

%   Define AIF estimation estting
opt.AIF.popAIF              = [];                       % container for population AIF (empty set will try to estimate AIF from data)
opt.AIF.AIFROI              = [];                       % keeping this empty will enforce ROI estimation
opt.AIF.AIFROItype          = 'vessel';                 % 'vessel' | 'wholebrain'
opt.AIF.lambda              = [];                       % regulariztion parameter to enforce vessel parameter, if non-empty vessel vp is allowed to follow a Gaussian distribution
opt.AIF.extractFrom         = 'abs';                    % 'abs' | 'real' | 'phase' | 'complex'
opt.AIF.integration         = 'sum';                    % 'sum' | 'trapz' | 'adaptive'
opt.AIF.paramterization     = [];                       % [] | '2Gammas+Sigmoid'
opt.AIF.initialParameters   = [];                       % []
opt.AIF.bloodT1             = [];                       % [] (will use measured T1) | some fixed value that is used
opt.AIF.bloodT2             = 0.275;                    % [] (will use measured T2 - currenlty not implemented) | some fixed value that is used
    % default values for blood T1 and T2star are given in Ref 3:
    %   T1  = 1932ms;
    %   T2* = 275ms;
opt.AIF.bloodr1             = 3.6;                      % Numbers taken from [4,5] for Gadovist
opt.AIF.bloodr2             = 6.3;
opt.AIF.bloodchi            = 320e-9;                   % molar susceptibility of the contrast agent: (roughly) 320 ppm/M = 320e-6 M^(-1) = 320e-9 L/mmol = 320e-9 1/mM [2]
opt.AIF.bloodHct            = 0.45;                     % average Hct value in US adult population according to Ref. [1],[6]
opt.AIF.shifts              = [0 1 -1];                 % time bin shifts that are applied prior to fitting, best fit is chosen
                                                        %   decimal shifts are Akima interpolated
                                                        %   postive shifts AIF right
                                                        %   negative shifts AIF left


%   Define Direct Fitting settings
opt.DIRECT  = [];

%   Define MOCCO settings
opt.MOCCO.lambda                = 1e3;
opt.MOCCO.enforceVesselParam    = true;
opt.MOCCO.boundVesselParam      = [];
opt.MOCCO.compression           = [];       % if non-empty this will use RSVD to save time on the SENSE part of the reconstruction

% opt.MOCCO.compression.rank            = 8;
% opt.MOCCO.compression.num_of_samples  = 12;   % put twice the rank or something
% opt.MOCCO.compression.sampleDist      = 'uniform';
% opt.MOCCO.compression.replacement     = false;

%   Define MOCCO SENSE options
opt.SENSE.solver            = 'pcg';
opt.SENSE.MaxIter           = 30;
opt.SENSE.tol               = 1e-4;
opt.SENSE.kspaceDataType    = 'list';   %   'full'|'list'

%   Define Motion Correction options
opt.MC.interframe           = false;

end