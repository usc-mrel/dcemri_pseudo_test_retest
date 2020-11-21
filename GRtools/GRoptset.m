function GRopt = GRoptset(flag, app)
%   Set options for dynamic reconstruction
%   
%   Author: RML
%   Date: 04/2014
%   
%   Usage: GRopt = GRoptset
%   
%   Input:
%   flag: bitmask for recon type:
%           1:  Static
%           2:  Pre/Post
%           4:  View-sharing
%           8:  SparseSENSE dynamic
%           16: Angio
%           32: Motion corection
%   app: application: 'perm', 'perf', 'ang', 'anat'
%   
%   Output:
%   GRopt: options structure

if nargin < 1 || isempty(flag)
    flag = 1;
    warning('Setting flags for anatomic recon only');
end
if nargin < 2 || isempty(app)
    app = 'perm';
    warning('Setting application to perm');
end

%   Set defaults
GRopt.DO_MOTION_CORR  = 0;
GRopt.DO_STATIC       = 0;      %   Should be on
GRopt.DO_PRE_POST     = 0;
GRopt.DO_VIEW_SHARING = 0;
GRopt.DO_SPSENSE      = 0;
GRopt.DO_ANGIO        = 0;
GRopt.DO_MOCCO        = 0;
GRopt.DO_DIRECT       = 0;

%   These are just placeholders to make the structure easier to read
GRopt.CC        = [];
GRopt.MC        = [];
GRopt.ST        = [];
GRopt.PP        = [];
GRopt.VS        = [];
GRopt.SP        = [];
GRopt.MOCCO     = [];
GRopt.DIRECT    = [];

%   Define which recons to perform
if bitget(flag,1)
    GRopt.DO_STATIC = 1;
    disp('Static recon:       ON  (01)');
else
    disp('Static recon:       OFF (01)');
end
if bitget(flag,2)
    GRopt.DO_PRE_POST = 1;
    disp('Pre/post recon:     ON  (02)');
else
    disp('Pre/post recon:     OFF (02)');
end
if bitget(flag,3)
    GRopt.DO_VIEW_SHARING = 1;
    disp('View-sharing recon: ON  (04)');
else
    disp('View-sharing recon: OFF (04)');
end
if bitget(flag,4)
    GRopt.DO_SPSENSE = 1;
    disp('SparseSENSE recon:  ON  (08)');
else
    disp('SparseSENSE recon:  OFF (08)');
end
if bitget(flag,5) && bitget(flag,4) && bitget(flag,2)
    GRopt.DO_ANGIO = 1;
    disp('Angio recon:      ON  (16)');
else
    disp('Angio recon:      OFF (16)');
end
if bitget(flag,6)
    GRopt.DO_MOTION_CORR = 1;
    disp('Motion correction:    ON  (32)');
else
    disp('Motion correction:    OFF (32)');
end
if bitget(flag,7)
    GRopt.DO_MOCCO = 1;
    disp('MOCCO reconstruction:  ON  (64)');
else
    disp('MOCCO reconstruction:  OFF (64)');
end


%%%%   Define application specific options    %%%%
switch lower(app)
    case {'t1'}
        %   Permeability
        disp('Setting parameters for T1 mapping');
        GRopt.SP.lambda             = 0.1;
        GRopt.SP.MaxIter            = 500;
        GRopt.SP.compression_levels = 25;
        GRopt.SP.segtype            = [];
        GRopt.SP.nbl                = 1;  %   Up to 24
        GRopt.SP.fast               = [];
        GRopt.SP.BF                 = [];
        
        GRopt.DIRECT.MaxIter        = 500;
        
        GRopt.MOCCO.lambda          = 5e4;
        GRopt.MOCCO.MaxIter         = 50;
        GRopt.MOCCO.senseTol        = 1e-6;
        
        GRopt.mtrxL = [];
        GRopt.mtrxH = [300 300 160];
        GRopt.binE = [];
        GRopt.app = 'T1Mo';
    case {'perm','permeability'}
        %   Permeability
        disp('Setting parameters for permeabiliy');
        GRopt.SP.lambda = [0 1e-3 0 0 1e-2 1e-2 0];
        GRopt.SP.lambda = [0 1 0 0 1 1 0];
            % these lambda values are initial lambda values and will be
            % re-adjusted to achieve the desired compression levels
            % specified in GRopt.SP.compression_levels
            % lambdas:
            % 1) GRAPPA
            % 2) complex Wavelet
            % 3) finte difference based TV
            % 4) Fourier domain filter
            % 5) View-sharing
            % 6) Proximity to reference image (experimental: can use the VS
            % image as reference)
            % 7) wavelet on difference to reference image
        GRopt.SP.compression_levels = [0 20 0 0 30 30 0];
        GRopt.SP.segtype = 1;
        GRopt.SP.nbl = 8;  %   Up to 24
        GRopt.SP.fast = 0;
        GRopt.SP.BF = [];
        GRopt.SP.MaxIter            = 200;
        
        GRopt.MOCCO.lambda          = 1e4;
        GRopt.MOCCO.MaxIter         = 100;
        
        GRopt.mtrxL = [];
        GRopt.mtrxH = [256 240 120];
        GRopt.binE  = 5; % set to 5s bins 
        GRopt.app = 'perm';
    case {'perf','perfusion'}
        %   Perfusion
        disp('Setting parameters for perfusion');
        GRopt.SP.lambda = [0 1e-3 0 0 1e-2 1e-2 0];
        GRopt.SP.compression_levels = [0 20 0 0 30 30 0];
        GRopt.SP.segtype = 1;
        GRopt.SP.nbl = 1;  %   Up to 16
        GRopt.SP.fast = 0;
        GRopt.SP.BF = [];
        
        GRopt.mtrxL = [];
        GRopt.mtrxH = [256 240 120];
        GRopt.binE = [];
        GRopt.app = 'perf';
    case {'ang','angio','angiography'}
        %   Angiography
        disp('Setting parameters for angiography');
        GRopt.SP.lambda = [0 1e-3 0 0 1e-3 1e-3 0];
        GRopt.SP.compression_levels = [0 20 0 0 30 30 0];
        GRopt.SP.segtype = 1;
        GRopt.SP.nbl = 1;  %   Up to 32
        GRopt.SP.fast = 0;
        GRopt.SP.BF = [];
        
        GRopt.mtrxL = [];
        GRopt.mtrxH = [256 240 120];
        GRopt.binE = [];
        GRopt.app = 'angio';
    case {'anat','anatomical'}
        %   Anatomical
        disp('Setting parameters for anatomical imaging');
        GRopt.SP.lambda = [0 1e-3 0 0 1e-2 1e-2 0];
        GRopt.SP.compression_levels = [0 40 0 0 75 75 0];
        GRopt.SP.segtype = 1;
        GRopt.SP.nbl = 1;
        GRopt.SP.fast = 0;
        GRopt.SP.BF = [];
        
        GRopt.mtrxL = [];
        GRopt.mtrxH = [256 240 120];
        GRopt.binE = [];
        GRopt.app = 'anat';
    otherwise
        error('Unknown application');
end
%%%%    END APPLICTION SPECIFIC OPTIONS    %%%


%%%%    Set common options    %%%%
GRopt.Dobj                  = 0.001;
GRopt.DCM                   = false;  %   DICOM output
GRopt.DCM_install           = false;
GRopt.parallel              = 1;
GRopt.class                 = 'single';
GRopt.bolus_arrival         = -1;     %   -1 for auto bolus detection
GRopt.variable_fr           = 0;      %   Flag for variable frame rate
GRopt.debug                 = [];     % select PE planes you want to use for debugging. The data will be subsampled respectively.
GRopt.gpuid                 = [];     % set to id of GPU if you want to pin it down
GRopt.JobStorageLocation    = [];     % folder for matlab to store the local files of the parallel pool

%   Define motion correction options
GRopt.MC.interframe = false;
GRopt.MC.mtrxL      = [48 48 48];
GRopt.MC.mtrxH      = [];
GRopt.MC.binE       = GRopt.binE;

%   Define static image options
GRopt.ST.lambda = [0 1e-2 0 0 0 0 0];
GRopt.ST.compression_levels = [0 30 0 0 0 0 0];

%   Define view sharing options
GRopt.VS.fwhm = 2.0;

%   Define pre/post contrast options
GRopt.PP.lambda = [0 1e-2 0 0 0 0 0];
GRopt.PP.compression_levels = [0 30 0 0 0 0 0];
GRopt.PP.windows = [180 180];

%   Define coil compression option
%   Set to [] for automatic
%   Set large (i.e., Inf) for none
GRopt.CC.nr = [];

% Set default VFA parameters:
GRopt.VFA.nFA = 6;                                              % number of flip angle prior to DCE acquisition flip angle
GRopt.VFA.nTR = 5000*GRopt.VFA.nFA;                             % number of TRs for VFA
GRopt.VFA.FAs = logspace(log10(1), log10(25), GRopt.VFA.nFA+1); % the current implementation uses log distributed flip angles

% Set names of files that need to be loaded
GRopt.filename.pre_post = 'GRrecon_perm_pre_post.mat';
GRopt.filename.t1map    = 'GRrecon_T1Mo.mat';                       % T1/M0 maps for DCE modelling
GRopt.filename.spsense  = ['GRrecon_' GRopt.app '_SPSENSE.mat'];    % 
