function [imgR, TK, OUT] = GRestimation(GRopt,path, outputfile_suffix)
%function [imgR, TK, OUT] = GRestimation(GRopt,path, outputfile_suffix)
%
%   Routine to estimate TK parameters
%
%   Inputs:
%       GRopt               structure for reconstruction/estimation options
%       path                absolut path to data, i.e., to p-file
%                           set to [], for estimation in current folder
%       outputfile_suffix   suffix to attach to output files (this might be
%                           useful if different parameter settings are
%                           supposed to be explored)
%
%   Outputs:
%       imgR                anatomic image time series
%       TK                  container of TK parameters
%       OUT                 additional output information
%                               .patient        patient info
%                               .bolus_arrival  bolus_arrival in seconds
%                                               (relative to begin of scan)
%                               .binE           vector of time bin edges
%                               .TKESTOUT       container of intermediate
%                                               auxiliary results of TK estimation
%
%   Usage:
%       GRopt = GRoptset(64, 'perm');
%       GRopt.mtrxH = [256 240 120];
%       GRopt.DO_MOCCO = 0;
%       GRopt.DO_SPSENSE = 1;
%       [imgR, TK, OUT] = GRestimation(GRopt, []);
%
%
%   Author
%       Yannick Bliesener, bliesene@usc.edu, MREL at USC, 2019/2020
%       R. Marc Lebel, mlebel@gmail.com
%


%   Start clock and open log file
stic = tic;
if nargin < 3 || isempty( outputfile_suffix )
   logfile = ['GRestimationlog_' GRopt.app '.txt'];
else
    logfile = ['GRestimationlog_' GRopt.app '_' outputfile_suffix '.txt'];
end
fid = fopen(logfile,'w');
fprintf(fid,'--- Cartesian radial image reconstruction ---\n');
fprintf(fid,'R. Marc Lebel\nmlebel@gmail.com\n\n');
fprintf(fid,'Yannick Bliesener\nbliesene@usc.edu\n\n');

%   Use current directory as default
if nargin<2 || isempty(path)
    path = cd;
end
cd(path);
fprintf(fid,'Date:                %s\n',datestr(now));
fprintf(fid,'Directory:           %s\n',path);

OUT = [];

matfile_version = '-v7.3';

% read current git commit hash
% this githash should be written to all data files that are produced so
% that the source code can be retrieved from the data file
githash = [];
repofolder = fileparts(mfilename('fullpath'));
[status,cmdout] = system(['(cd ' repofolder '; git rev-parse --verify HEAD)']);
if ~status
    githash = strtrim(string(cmdout));
else
    fprintf(fid,'Could not read GIT commit hash!\n');
end
clear status cmdout


%   Read raw data
fprintf(fid,'\n\n---- PREP ----\nReading raw data...');tic;
[k,petab,hdr,tstamp,pfle] = GRread2(path,GRopt.mtrxL,GRopt.mtrxH);
[np, nv, ns, nr] = get_key_params2(pfle);
tr = hdr.ImageData.tr / 1e6;
te = hdr.ImageData.te / 1e6;
voxel_size = [hdr.ImageData.dfov/ np, hdr.ImageData.dfov / nv, hdr.ImageData.slthick];
FOV = [hdr.RawHeader.fov, hdr.RawHeader.fov*hdr.RawHeader.phase_scale, hdr.ImageData.slthick*ns];
OUT.patient.weight = hdr.ExamData.patweight;
OUT.patient.sex = hdr.ExamData.patsex;
OUT.patient.age = hdr.ExamData.patage;
switch GRopt.class
    case 'double'
        k = double(k);
    case 'single'
        k = single(k);
end
fprintf(fid,'done (%g sec)\n',toc);


%   Perform coil compression
if isempty(GRopt.CC.nr) || GRopt.CC.nr < nr
    fprintf(fid,'Performing coil compression...');tic;
    [k, nr] = geometric_coil_compression(k,GRopt.CC.nr,petab);
    fprintf(fid,'done (%g sec)\n',toc);
    fprintf(fid,'\tCompressed to %g virtual coils\n',nr);
    GRopt.CC.nr = nr;
    pfle.channels = nr;
end

% Read flip angle table (Yannick)
dirs = dir('*fliptable.txt');
fliptable = [];
if length(dirs) == 1
    fprintf(fid,'Reading flip angle table ...');tic;
    fliptable = dlmread( dirs.name );
        % load the flip angle table
    FAs = unique(fliptable);
        % the FAs for VFA
    nFA = length(FAs) - 1;
        % number of flip angle prior to DCE acquisition flip angle
    nTR = find( (fliptable == fliptable(end)), 1, 'first') + 500; % just put a buffer of 1000 TRs for the last flip angle
        % number of TRs for VFA, since the data gets cut later this remains
        % really just an estimate and technically needs to be
        % updated ... but it's good enough to offset the bolus
        % arrival detection!
    fprintf(fid,'done (%g sec)\n',toc);
    
    % just overwrite the defaul values
    GRopt.VFA.nFA = nFA;
    GRopt.VFA.nTR = nTR;
    GRopt.VFA.FAs = FAs;
else
    warning('Using default flip angle table!')
    fprintf(fid,'Could not find flip angle table!\n')
    fprintf(fid,'->Using default flip angle table!\n')
    
    nFA = GRopt.VFA.nFA;
        % number of flip angle prior to DCE acquisition flip angle
    nTR = GRopt.VFA.nTR;
        % number of TRs for VFA
    FAs = GRopt.VFA.FAs;
     % the current implementation uses log distributed flip angles
end


%   Remove dummy pulses
% NOTE: STARDCE probably won't have that
if hdr.ImageData.user46
    fprintf(fid,'Removing dummy pulses...\n');tic;
    NSTST = hdr.ImageData.user46;
    k = k(:,(NSTST+1):end,:,:,:);
    petab = petab((NSTST+1):end,:);
    tstamp = tstamp((NSTST+1):end);
    if ~isempty(fliptable)
        fliptable = fliptable((NSTST+1):end);
    end
    fprintf(fid,'done (%g sec)\n',toc);
else
    NSTST = 0;
end

%  Remove Eddy current calibration region
% NOTE: STARDCE probably won't have that
if hdr.ImageData.user47
    fprintf(fid,'Performing eddy current correction...\n');
    fprintf(fid,'\tExtracting data...');tic;
    NEDDY = hdr.ImageData.user47;
    kE = k(:,1:NEDDY,:,:,:);
    k = k(:,(NEDDY+1):end,:,:,:);
    petab = petab((NEDDY+1):end,:);
    tstamp = tstamp((NEDDY+1):end);
    if ~isempty(fliptable)
        fliptable = fliptable((NEDDY+1):end);
    end
    fprintf(fid,'done (%g sec)\n',toc);
else
    NEDDY = 0;
end
 
%   Perform eddy current estimation
%   TODO: should adjust which receivers to use if CC not used
% fprintf(fid,'\tEstimating correction terms...');tic;
% OUT.EC = GReddy2(kE,1);
% fprintf(fid,'done (%g sec)\n',toc);
% fprintf(fid,'\tApplying correction terms...');tic;
% k = GReddy_apply(OUT.EC,k,petab);
% fprintf(fid,'done (%g sec)\n',toc);

%   Extract sensitivity region
if hdr.ImageData.user48
    % NOTE: STARDCE probably won't have that
    NCAL = hdr.ImageData.user48^2;
    kS = GRreshape2(k,pfle,petab,[(NCAL+NEDDY+NSTST)*tr-tr, 1e6],tstamp,0);
    kS = kS(:,:,:,1,:);
else
    % Condensing (averaging) all temoral data into a single frame to get a
    % sensitivity map
    kS = GRreshape2(k,pfle,petab,Inf,tstamp,0);
    NCAL = 0;
end

%   Strip T1 mapping
% stripping VFA, coil sensitivty calibration region, eddy current
% calibration region, and dummy pulses
indstart = find(tstamp>(nTR + NCAL + NEDDY + NSTST)*tr,1);
fprintf(fid,'Removing T1 mapping data (0 - %3.2f sec)\n',tstamp(indstart));
k = k(:,indstart:end,:,:,:);
petab = petab(indstart:end,:);
tstamp = tstamp(indstart:end);
if ~isempty(fliptable)
    fliptable = fliptable(indstart:end);
end
nviews = size(k,2);

%   Discard some time at the beginning of each pass
tss = 2.00;
ind = zeros(size(tstamp));
passdur = tstamp(end)/pfle.passes;
for i = 1:pfle.passes
    ind = ind + (tstamp>(i*passdur) & tstamp<=(i*passdur+tss));
end
ind = ~ind;
k = k(:,ind,:,:,:);
tstamp = tstamp(ind);
petab = petab(ind,:);
if ~isempty(fliptable)
    fliptable = fliptable(ind);
end
clear ind tss i


%   Detect bolus arrival
if GRopt.bolus_arrival < 0
    fprintf(fid,'Determining bolus arrival...');tic;
    indcent = (petab(:,1)==nv/2+1 & petab(:,2)==ns/2+1);
        % VFA part has been stripped
    %kcent = k(np/2+1,indcent,:,:,:);
    kcent = k(:,indcent,:,:,:);
    kcent = fftshift(iFastFT(kcent,1,1),1);
    kcent = sqrt(sum(abs(kcent).^2,5));
    tcent = tstamp(indcent);
    OUT.bolus_arrival = GRestbolus(kcent,tcent);
    clear indcent kcent tcent;
    fprintf(fid,'done (%g sec)\n',toc);
    if isempty(OUT.bolus_arrival)
        fprintf(fid,' Bolus not detected...aborting\n');
        fclose(fid);
        error('Bolus not detected');
    else
        fprintf(fid,' Bolus detected at %05.1f s\n',OUT.bolus_arrival);
    end
else
    OUT.bolus_arrival = GRopt.bolus_arrival;
    if isempty(OUT.bolus_arrival) || OUT.bolus_arrival == 0 || OUT.bolus_arrival > tstamp(end)
        fprintf(fid,' Unknown/no bolus arrival\n');
        fclose(fid);
        error('Unknown/no bolus arrival');
    else
        fprintf(fid,' Bolus specified at %05.1f s\n',OUT.bolus_arrival);
    end
end
if ~isempty(OUT.bolus_arrival)
    GRopt.PP.windows(1) = OUT.bolus_arrival;
    GRopt.PP.windows(2) = OUT.bolus_arrival;
end


%   Set frame rate
OUT.binE = GRsetbinE(OUT.bolus_arrival,GRopt,tstamp);
fprintf(fid,'Temporal bins (s):\n');
fprintf(fid,' %05.1f',OUT.binE);
fprintf(fid,'\n');


%   Estimate coil sensitvities
fprintf(fid,'Estimating coil sensitivities...');tic
[S,sf] = GRcoilsens(kS,NCAL);
k = sf*k;
clear kS
fprintf(fid,'done (%g sec)\n',toc);

% Load & Clean T1/M0 maps
load(GRopt.filename.t1map, 'B1', 'T1', 'Mo', 'imgFA', 'sf2');
k = sf2 * k;
sf2 = 1;    % important to prevent it from being reapplied during phase correction

% load debugData.mat
GRopt.fid = fid;

% determine brain mask
load(GRopt.filename.pre_post, 'imgPOST')

if 0
%     refImg = imgFA(:,:,:,end);
    refImg = imgPOST;
    % pick last flip angle because this one is sharpest
    
    brainmask = abs(refImg) > 1e4;
        % just to generate an initial brainmask that is overwritten
    for i=1:size(brainmask,1)
        tmp = refImg(i,:,:);
        maxSig = prctile(abs(tmp(:)), 75);
        T = graythresh(min(1, abs(refImg) ./ maxSig) ) * maxSig;
        brainmask(i,:,:) = abs(tmp) > T;
    end
%     clear imgFA

    % detect slices without any signal
    tmp = refImg(:);
    maxSig = prctile(abs(tmp(:)), 75);
    T = graythresh(abs(refImg) ./ maxSig ) * maxSig;
    tmp = (refImg < T);
    tmp = sum(sum(tmp,2),3);
    
    ind = (tmp == (size(brainmask,2)*size(brainmask,3)));
    brainmask(ind, :, :) = false;

    clear tmp ind T maxSig
else
%     maxM0 = prctile(abs(Mo(:)), 75);
%     T = graythresh(abs(Mo) ./ maxM0) * maxM0;
%     brainmask = abs(Mo) > T;
%     
%     imgFA = imgFA(:,:,:,end);
    refImg = imgPOST;
    
    maxSig = prctile(abs(refImg(:)), 75);
    T = graythresh(min(1, abs(refImg) ./ maxSig)) * maxSig;
    brainmask = abs(refImg) > T;
end
% brainmask = bwareaopen(brainmask, 9);
for i=1:size(brainmask,1)
%     brainmask(i,:,:) = imfill(squeeze(brainmask(i,:,:)), 'holes');
    brainmask(i,:,:) = imdilate(squeeze(brainmask(i,:,:)), strel('disk', 1));
end

% clean M0
variableInfo = who('-file', GRopt.filename.t1map);
    % check if preprocessed T1 maps have already been stored to file for
    % speed reasons
if ~ismember('R1', variableInfo)
    Mo( ~brainmask(:) ) = 0;
    R1 = 1./T1 * 1000;

    ind = T1 > 1e4;
    R1(ind) = 0;
    Mo(ind) = 0;
    % brainmask(ind) = false;
    for i = find(ind(:))'
        [x,y,z] = ind2sub(GRopt.mtrxH, i);
        if (x==1) || (x == GRopt.mtrxH(1)) || (y == 1) || (y==GRopt.mtrxH(2)) || (z == 1) || (z==GRopt.mtrxH(3))
            continue
        end
        [X, Y, Z] = meshgrid(x-1:x+1, y-1:y+1, z-1:z+1);

        tmp = R1(X(:), Y(:), Z(:));
        R1(i) = median(tmp(:));

        tmp = Mo(X(:), Y(:), Z(:));
        Mo(i) = median(tmp(:));
    end

    ind = isinf(R1);
    R1(ind) = 0;
    Mo(ind) = 0;

    ind = isinf(Mo);
    R1(ind) = 0;
    Mo(ind) = 0;

    clear T1
else
    % load preprocessed R1 and Mo
    fprintf(fid,'Loading preprocessed R1 maps from file.\n');
    load(GRopt.filename.t1map, 'R1', 'Mo');
end


% for debugging only: downsample
if ~isempty(GRopt.debug)
    ind = GRopt.debug;

    tmp_opt.FTshift = 1;
    tmp_opt.FTdim = 1;
    k = iFastFT(fftshiftF(k, 1),tmp_opt);
        % careful with the shift dimension
        % this is otherwise called in GRreshape2
    k = k(ind, :, :, :, :);
    k = fFastFT(k,tmp_opt);
    k = ifftshiftF(k,1);

    S = S(ind, :, :, :, :);
    R1 = R1(ind, :, :, :, :);
    Mo = Mo(ind, :, :, :);
    B1 = B1(ind, :, :, :);
    imgFA = imgFA(ind, :, :, :);
    brainmask = brainmask(ind, :, :);
    
    np = length(GRopt.debug);
    pfle.xRes = length(GRopt.debug);
    GRopt.mtrxH(1) = length(GRopt.debug);
    clear ind tmp_opt
    
    save('debugData.mat', '-regexp', '^(?!(fid)$).', '-v7.3')
end


%   Get generic Tk estimation recon options
opt              = TKESTIMATION_optset;
opt.Dobj         = GRopt.Dobj;
opt.class        = GRopt.class;
opt.FOV          = FOV;
opt.gpuid        = GRopt.gpuid;

if isfield(GRopt, 'MOCCO')
    for fn = fieldnames(GRopt.MOCCO)'
        opt.MOCCO.(fn{1}) = GRopt.MOCCO.(fn{1});
    end
end
if isfield(GRopt, 'TK')
    for fn = fieldnames(GRopt.TK)'
        opt.TK.(fn{1}) = GRopt.TK.(fn{1});
    end
end
if isfield(GRopt, 'AIF')
    for fn = fieldnames(GRopt.AIF)'
        opt.AIF.(fn{1}) = GRopt.AIF.(fn{1});
    end
end
if isfield(GRopt, 'SENSE')
    for fn = fieldnames(GRopt.SENSE)'
        opt.SENSE.(fn{1}) = GRopt.SENSE.(fn{1});
    end
end
if isfield(GRopt, 'MC')
    for fn = fieldnames(GRopt.MC)'
        opt.MC.(fn{1}) = GRopt.MC.(fn{1});
    end
end

%   General placeholder
TK = [];

% create output file
if nargin < 3 || isempty(outputfile_suffix)
    matname = ['GRrecon_TK_' GRopt.app '.mat'];
else
    matname = ['GRrecon_TK_' GRopt.app '_' outputfile_suffix '.mat'];
end
save(matname,'OUT', 'opt', 'githash', matfile_version);

if nargin < 3 || isempty(outputfile_suffix)
    matnameimg = ['GRrecon_' GRopt.app '.mat'];
else
    matnameimg = ['GRrecon_' GRopt.app '_' outputfile_suffix '.mat'];
end

%   Open parallel pool
ncores = feature('numCores');
parpl = gcp('nocreate');
pc = [];
if isempty(parpl)
    % create a local cluster object
    pc = parcluster('local');

    % explicitly set the JobStorageLocation
    % prevents race conditions on clusters
    % should go unnoticed if no slurm is used
    if getenv('SLURM_JOB_ID')
        slurmIDs = {getenv('SLURM_JOB_ID'), getenv('SLURM_ARRAY_TASK_ID')};
        slurmIDs = slurmIDs(~cellfun('isempty', slurmIDs));
        if isempty(GRopt.JobStorageLocation)
            job_folder = fullfile(pc.JobStorageLocation, strjoin(cat(2, {'job'}, slurmIDs),'_'));
        else
            job_folder = fullfile(GRopt.JobStorageLocation, strjoin(cat(2, {'job'}, slurmIDs),'_'));
        end
        evalc(['system(''mkdir -p ' job_folder ''')']);
        pc.JobStorageLocation = job_folder;
    end        
    
    parpl = parpool(pc, min([np ncores]), 'IdleTimeout', Inf);
    
end


%% %%%% MOCCO TK ESTIMATION %%%%%%
imgR = [];
if GRopt.DO_MOCCO
    fprintf(fid,'\n\n---- MOCCO TK ESTIMATION ----\n');
    opt.S = S;
    opt.M0 = Mo;
    opt.B1 = B1;
    opt.R1 = R1;
    opt.FA = FAs(end);
    opt.TR = tr;
    opt.TE = te;
    opt.voxel_size = voxel_size;
    opt.phase_corr = [1 2 3];
    opt.Ph_cf = [1 1 1];
    opt.brainmask = brainmask;
    
    %   Recon image
    [imgR, TK, TKESTOUT] = MTKrecon(k,opt,GRopt,petab,pfle,OUT.binE,OUT.bolus_arrival,tstamp);
    OUT.MOCCO = TKESTOUT;
    
    %   Save
%     fprintf(fid,'Phase correcting...');tic;
%     imgR = PHcorrection(imgR,sf2);
%     %imgR = flip(imgR,3);
%     fprintf(fid,'done (%g sec)\n',toc);
    fprintf(fid,'Saving .mat...');tic;
    save(matname, 'TK', 'TKESTOUT', '-append');
    save(matnameimg, 'imgR', 'githash', matfile_version);
    fprintf(fid,'done (%g sec)\n',toc);
    if GRopt.DCM
        fprintf(fid,'Writing DICOMs...');tic;
        GRwrite_dicom(imgR,pfle,hdr,'ImgDY','imgDY','MOCCO dynamic',4,GRopt.DCM_install);
        fprintf(fid,'done (%g sec)\n',toc);
    end
%     clear imgR
end


%% %%%% TK ESTIMATION %%%%%%
if GRopt.DO_SPSENSE
    fprintf(fid,'\n\n---- TK ESTIMATION ----\n');
%     opt.S = S;
    opt.M0 = Mo;
    opt.B1 = B1;
    opt.R1 = R1;
    opt.FA = FAs(end);
    opt.TR = tr;
    opt.TE = te;
    opt.voxel_size = voxel_size;
    opt.phase_corr = [1 2 3];
    opt.Ph_cf = [1 1 1];
    opt.brainmask = brainmask;
    
    % load image data
    load(GRopt.filename.spsense, 'imgR')
    load(GRopt.filename.t1map, 'sf2');
    imgR = sf2 .* imgR;
    
    %   Recon image
    [TK, TKESTOUT] = TKestimation(imgR,opt,GRopt,OUT.binE,OUT.bolus_arrival);
    OUT.TKESTOUT = TKESTOUT;
    
    %   Save
%     fprintf(fid,'Phase correcting...');tic;
%     imgR = PHcorrection(imgR,sf2);
%     %imgR = flip(imgR,3);
%     fprintf(fid,'done (%g sec)\n',toc);
    fprintf(fid,'Saving .mat...');tic;
    save(matname, 'TK', 'TKESTOUT', '-append');
    fprintf(fid,'done (%g sec)\n',toc);
end

%   Close parallel pool
delete(parpl);
if ~isempty(pc)
    pc.Jobs.delete;
end

%   End
save(matname,'OUT','GRopt','-append');
fprintf(fid,'\n ---- DONE ----\nReconstruction time: %g min\n\n',...
        0.01*round(100*toc(stic)/60));
fclose(fid);

return

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NESTED FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [imgR, TK, outputInfo] = MTKrecon(k,opt,GRopt,pe_in,pfle_in,binE_in,BAT_in,tstamp_in)
        
        %   Get key parameters
        [np_dy, nv_dy, ns_dy, nr_dy] = get_key_params2(pfle_in);
        
        fprintf(GRopt.fid,'Reshaping raw data...');tic;
        
        %   Set some options and clean up others to remove large unused variables
        opt.fid = GRopt.fid;
        opt.time = binE_in - BAT_in;
        [~, opt.ibolus] = min(abs(opt.time));
        opt.tbolus = BAT_in;
        
        if isempty(opt.MOCCO.compression)
            
            %   Rebin data and shift
            [kdy,Udy] = GRreshape2(k,pfle_in,pe_in,binE_in,tstamp_in,0);
            ph = size(kdy,4);
            R = nv_dy*ns_dy/median(sum(sum(squeeze(Udy(1,:,:,:,1)),1),2));
            
            for iind = 1:nr_dy
                kdy(:,:,:,:,iind)  = iFastFT(kdy(:,:,:,:,iind),1,1);
                kdy(:,:,:,:,iind)  = fftshift(fftshift(kdy(:,:,:,:,iind),2),3);
                Udy(:,:,:,:,iind)  = fftshift(fftshift(Udy(:,:,:,:,iind),2),3);
            end
            clear iind jind
            kdy = kdy(Udy);
            kdy = reshape(kdy,[np_dy numel(kdy)/np_dy]);
            
            opt.size = [np_dy nv_dy ns_dy ph nr_dy];
            opt.U = find(Udy(:,:,:,:,:));
            
        else
            % in the long run this is section is the way to go
            % it performs the same operations as the other section yet does
            % it without allocating the full k-space across all time frames
            % hence it scales better in memory efficiency when more time
            % frames are used
            
            %   Rebin data and shift
            [kdy,opt.U, ph] = GRprepkspace(k,pfle_in,pe_in,binE_in,tstamp_in,0);
            
            opt.size = [np_dy nv_dy ns_dy ph nr_dy];
        end

        fprintf(GRopt.fid,'done (%g sec)\n',toc);
        fprintf(GRopt.fid,'Lambda = [');
        fprintf(GRopt.fid,' %g',opt.MOCCO.lambda);
        fprintf(GRopt.fid,']\n');
        
        [TK,imgR,outputInfo] = MOCCOTK(kdy,opt);
        
    end

    function [TK, outputInfo] = TKestimation(img,opt,GRopt,binE_in,BAT_in)
        
        outputInfo = struct();
        
        %   Set some options and clean up others to remove large unused variables
        opt.fid = GRopt.fid;
        opt.time = binE_in - BAT_in;
        [~, opt.ibolus] = min(abs(opt.time));
        opt.tbolus = BAT_in;
        opt.size = size(img);
        
        if ~isempty(GRopt.gpuid) && (GRopt.gpuid > 0)
            gpuDevice(GRopt.gpuid);
        end
        
        % Re-estimate phase correction term
        opt.Ph = repmat(exp(sqrt(-1)*angle(img(:,:,:,1))), [1 1 1 opt.size(4)]);
        img = opt.iPH(img,opt);
        
        % estiamte AIF
        fprintf(GRopt.fid,'Estimating AIF...');tic;
        [AIF, AIFROI, AIFparam, AIFtheta, vessels] = estimateAIF(img, opt);
        fprintf(GRopt.fid,'done (%g sec)\n',toc);
        
        
        if isempty(AIF) || any(isnan(AIF)) || any(isinf(AIF))
            TK = [];
            outputInfo.AIFs = AIF;
            outputInfo.AIFROI = AIFROI;
            outputInfo.AIFparam = AIFparam;
            outputInfo.AIFtheta = AIFtheta;
            outputInfo.vessels = vessels;
            fprintf(GRopt.fid,'\tFailed to estimate AIF\n');
            return
        end        
        
        img = abs(img);
        
        % Motion correction
        tform = [];
        if opt.MC.interframe
            fprintf(GRopt.fid,'Motion correction...');tic;
            [img,tform,tformI,Mot] = registervolumes(img ,opt.FOV);
            fprintf(GRopt.fid,'done (%g sec)\n',toc);
        end
        
        % Signal to R1 map
        R1 = iSPGR(img, opt);
        
        % R1 to CTC
        conc = iFXL(R1, opt);
        conc( isinf(conc) | isnan(conc) ) = 0;
    
        %% MOCCO: TK estimation
        fprintf(GRopt.fid,'Estimating TK parameters...');tic;
        initial_param = [];
        [TK, ~, FITOUT] = estimateTK(conc, AIF, initial_param, opt);
        outputInfo = FITOUT;
        fprintf(GRopt.fid,'done (%g sec)\n',toc);
        
        outputInfo.AIFs = AIF;
        outputInfo.AIFROI = AIFROI;
        outputInfo.AIFparam = AIFparam;
        outputInfo.AIFtheta = AIFtheta;
        outputInfo.vessels = vessels;
        outputInfo.tform = tform;
        
    end


end