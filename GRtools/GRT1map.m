function [T1, Mo, OUT] = GRT1map(GRopt, path, outputfile_suffix)
%function [T1, Mo, OUT] = GRT1map(GRopt, path, outputfile_suffix)
%
%   Function to estimate baseline T1/M0 maps from VFA data
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
%       T1                  baseline T1 map
%       Mo                  proton density map
%       OUT                 additional output information
%
%   Usage:
%       GRopt = GRoptset(8, 't1');
%       GRopt.mtrxH = [256 240 120];
%       [T1, Mo, OUT] = GRT1map(GRopt,[]);
%
%
%   Author
%       R. Marc Lebel, mlebel@gmail.com
%       Zhibo Zhu, zhibozhu@usc.edu
%       Yannick Bliesener, bliesene@usc.edu, MREL at USC, 2019/2020
%       
%

%   Start clock and open log file
stic = tic;
if nargin < 3 || isempty( outputfile_suffix )
    logfile = ['GRT1log_' GRopt.app '.txt'];
else
    logfile = ['GRT1log_' GRopt.app '_' outputfile_suffix '.txt'];
end
fid = fopen(logfile,'w');
fprintf(fid,'--- Cartesian radial image reconstruction ---\n');
fprintf(fid,'R. Marc Lebel\nmlebel@gmail.com\n\n');
fprintf(fid,'Z. Zhu \nzhibozhu@usc.edu\n\n');
fprintf(fid,'Y. Bliesener \nbliesene@usc.edu\n\n');

%   Use current directory as default
if nargin<2 || isempty(path)
    path = cd;
end
cd(path);
fprintf(fid,'Date:                %s\n',datestr(now));
fprintf(fid,'Directory:           %s\n',path);

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

%   Get B1 map
fprintf(fid,'\n\n---- PREP ----\n');
fprintf(fid,'Getting B1+ map\n');
dirs = dir('../B1*');
if length(dirs) == 1
    pname{1} = ['../' dirs(1).name];
else
    pname{1} = uigetdir(pwd,'Select B1+ dicom folder');
end
if pname{1} ~= 0
    [B1, B1param] = getB1map(pname);
    switch GRopt.class
        case 'double'
            B1 = double(B1);
        case 'single'
            B1 = single(B1);
    end
    fprintf(fid,'\tB1+ maps read from %s\n',pname{1});
    B1rot = [B1param.ImageOrientationPatient(1:3,1), B1param.ImageOrientationPatient(4:6,1)];
    % these are the first two columns of the rotation matrix to get from
    % display coordinates (r') to physical coordinates (r: x = R/L; y = A/P; z = S/I)
    % NOTE: all dicom images of the B1+ map have to orientated identically
    % r = B1rot * r'
    B1rot = [B1rot, zeros(3,1)];
    if B1rot(1,1) == 1 && B1rot(2,2) == 1
        % Keck STARDCE setup
        B1rot(3,3) = 1;
    elseif B1rot(1,1) == 1 && B1rot(3,2) == -1
        % Calgary setup
        B1rot(2,3) = 1;
    else
        fprintf(fid,'\tUnkown orientation of B1 map\n');
        error('Unkown orientation of B1 map')
    end
else
    B1 = ones([2 2 2],GRopt.class);
    fprintf(fid,'\tNo B1 map selected\n');
    B1rot = [];
end
clear F nsB1 dirs pname B1param

%   HACK. Not sure why this is needed. T1 values come out way off otherwise.
%   Need to figure this out.
%B1 = 1.25*B1;

%   Read raw data
fprintf(fid,'Reading raw data...');tic;
[k,petab,hdr,tstamp,pfle] = GRread2(path,GRopt.mtrxL,GRopt.mtrxH);
[np, nv, ns, nr] = get_key_params2(pfle);
tr = hdr.ImageData.tr / 1e6;
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
    fprintf(fid,'Could not find flip angle table!\n');
    fprintf(fid,'->Using default flip angle table!\n');
    
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
    if ~isempty(fliptable)
        fliptable = fliptable((NSTST+1):end);
    end
    tstamp = tstamp((NSTST+1):end);
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
    if ~isempty(fliptable)
        fliptable = fliptable((NEDDY+1):end);
    end
    tstamp = tstamp((NEDDY+1):end);
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
    NCAL = [];
end

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
    fliptable = fliptable(ind, :);
end
clear ind tss i

% plot k-space center
indcent = (petab(:,1)==nv/2+1 & petab(:,2)==ns/2+1);
tcent = tstamp(indcent);
kcent = k(:,indcent,:,:,:);
kcent = sqrt(sum(abs(kcent).^2,5));
kcent = kcent(np/2+1,:);
kcent = medfilt1(kcent,3);
opt.kcent = kcent;
save kcent kcent tcent


%   Detect bolus arrival
if GRopt.bolus_arrival < 0
    fprintf(fid,'Determining bolus arrival...');tic;
    %indcent = (petab(:,1)==nv/2+1 & petab(:,2)==ns/2+1 & ...
    %           tstamp(:) > (hdr.ImageData.user30+NCAL+NEDDY+NSTST)*tr);
    indcent = (petab(:,1)==nv/2+1 & petab(:,2)==ns/2+1 & tstamp(:) > nTR*tr);
        % get center of k-space space after initial 6 VFAs for T1
        % estimation which each take 2000 TRs
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


% Set parameters for VFA method: 

% Set bin edges
if isempty(fliptable)
    %OUT.binE = nTR/nFA * tr * (1:nFA) + (NCAL+NEDDY)*tr;   %   The first bin has a calibration zone, also time shifted by eddy current cal
    OUT.binE = nTR/nFA * tr * (1:nFA) + (NCAL+NEDDY+NSTST)*tr;   %   The first bin has a calibration zone, also time shifted by extra pulses.
else
    FAbins = false(size(fliptable));
    FAbins(2:end) = abs( fliptable(2:end) - fliptable(1:end-1) ) > 0;
    OUT.binE = find(FAbins)' * tr;
    clear FAbins    
end
OUT.binE = [OUT.binE OUT.bolus_arrival 1e6];    %   Add the time between the last flip angle and bolus arrival
    % Yannick: binE stores bin edges. second to last is 15 degree prior to
    % contrast injection, last bin is after bolus arrival to end (1e6)
%OUT.FA = exparray(1.5,hdr.ImageData.mr_flip,nFA+1);
% OUT.FA = linspace(1.5,15,10);
OUT.FA = FAs;
    
fprintf(fid,'Temporal bins (s):\n');
fprintf(fid,' %05.1f',OUT.binE);
fprintf(fid,'\n');

% create output file
if nargin < 3 || isempty(outputfile_suffix)
    matname = 'GRrecon_T1Mo.mat';
else
    matname = ['GRrecon_T1Mo_' outputfile_suffix '.mat'];
end
save(matname,'OUT', 'githash', matfile_version);

%   Discard some time after each flip change
%tss = 4.00;
tss = 1.00;
ind = zeros(size(tstamp));
for i = 1:(length(OUT.binE)-2)
    ind = ind + (tstamp>=OUT.binE(i) & tstamp<=(OUT.binE(i)+tss));
end
ind = ~ind;
tstamp = tstamp(ind);
k = k(:,ind,:,:,:);
petab = petab(ind,:);
if ~isempty(fliptable)
    fliptable = fliptable(ind,:);
end
clear ind tss i


%   Sort data and trim down to only flip angle maps
fprintf(fid,'Sorting data...');tic;
[k,U] = GRreshape2(k,pfle,petab,OUT.binE,tstamp,0);
k = k(:,:,:,1:(nFA+1),:);
U = U(:,:,:,1:(nFA+1),:);
fprintf(fid,'done (%g sec)\n',toc);


%   Flip data in slice direction
% k  = flip(k,3);
% U  = flip(U,3);
% kS = flip(kS,3);


%   Estimate coil sensitvities
fprintf(fid,'Estimating coil sensitivities...');tic
[S,~] = GRcoilsens(kS,NCAL);
clear kS;
fprintf(fid,'done (%g sec)\n',toc);

% Interpolate B1 map
if ~isempty(B1rot)
    % change orientation of B1 map if necessary
    % Matlab has rigid3d rotation as of Matlab 2020 but I am too lazy to
    % upgrade so someone else can do that if they want to... in the
    % meantime this hack might do
    
    refR = [1 0 0; 0 0 1; 0 -1 0];
    if norm(refR(:) - B1rot(:)) < eps
        % Calgary setup
        % do nothing
    elseif norm(reshape(eye(3), [],1) - B1rot(:)) < eps
        % Keck setup -> make it identical to the Calgary setup, i.e.,
        % specifiy B1 as a stack of coronal slices
        B1 = permute(B1, [3 2 1]);
        B1 = flip(B1, 1);
    else
        fprintf(fid,'\tUn-interpretable orientation of B1 map - this requires additional implementation effort!\n');
        error('Un-interpretable orientation of B1 map - this requires additional implementation effort!')
    end
end
B1 = imresize3d(B1,[np nv ns]);
save(matname,'B1','-append');

% for debugging only: downsample
if ~isempty(GRopt.debug)
    ind = GRopt.debug;
    GRopt.SP.nbl = []; % Skip rebin when in debug mode, Zhibo 03/17/2020

    tmp_opt.FTshift = 1;
    tmp_opt.FTdim = 1;
    k = iFastFT(k,tmp_opt);
    k = k(ind, :, :, :, :);
    k = fFastFT(k,tmp_opt);
    
    S = S(ind, :, :, :, :);
    U = U(ind, :, :, :, :);
    B1 = B1(ind, :, :);
    opt.kcent = opt.kcent .* (length(ind) / np); 
    [np, nv, ns, ~, nr] = size(k);
    clear ind tmp_opt
    
    save('debugT1est', '-regexp', '^(?!(fid)$).', '-v7.3')
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
    
    parpl = parpool(pc, min([np GRopt.SP.nbl ncores]), 'IdleTimeout', Inf);
    
end

% prepare GPU device
if ~isempty(GRopt.gpuid) && (GRopt.gpuid > 0)
   opt.gpuid = gpuDevice(GRopt.gpuid);
end


% reconstruct
if GRopt.DO_SPSENSE
    
    %   Create options structure
    opt.class = 'single';
    opt.FA = OUT.FA;
    opt.tr = tr;
    opt.S = S;
    opt.B1 = B1;
    opt.phase_corr = [1 2 3];
    opt.Ph_cf = [1 1 1];
    opt.size = [np nv ns nFA+1 nr];
    opt.FTshift = 1;
    opt.FTdim = [1 2 3];
    opt.fid = fid;
    opt.FOV = [hdr.RawHeader.fov ...
        hdr.RawHeader.fov*hdr.RawHeader.phase_scale ...
        hdr.ImageData.slthick*ns];
    
    %   Prep data
    opt.U = find(U);clear U
    k = k(opt.U);
    
    % copy solver specific parameters
    opt.MaxIter = GRopt.SP.MaxIter;
    opt.lambda  = GRopt.SP.lambda;

    %   Call iterative recon
    fprintf(fid,'Performing iterative T1 mapping');tm = tic;
    %[T1,Mo,imgFA] = SPDESPOT1(k,opt);
    [T1,Mo,imgFA] = SPDESPOT1_temp(k,opt);
    fprintf(fid,'done (%g sec)\n',toc(tm));
    
else
    fprintf(fid,'ERROR: No T1 estimation procedure specified!');tic;
    fclose(fid);
    T1 = [];
    Mo = [];
    return
end


%   Close parallel pool
delete(parpl);
if ~isempty(pc)
    pc.Jobs.delete;
end

%   Scale Mo and T1
sf2 = 24576/prctile(Mo(:),99);
Mo = sf2 * Mo;
T1 = T1*1000;   %   In ms

%   Save
fprintf(fid,'Saving .mat...');tic;
save(matname,'T1','Mo','sf2','imgFA','-append');
fprintf(fid,'done (%g sec)\n',toc);
if GRopt.DCM
    fprintf(fid,'Writing DICOMs...');tic;
    GRwrite_dicom(Mo,pfle,hdr,'Mo','Mo','Mo',6,GRopt.DCM_install);
    GRwrite_dicom(T1,pfle,hdr,'T1','T1','T1',7,GRopt.DCM_install);
    fprintf(fid,'done (%g sec)\n',toc);
end

%   End
save(matname,'OUT','GRopt','-append');
fprintf(fid,'\n ---- DONE ----\nReconstruction time: %g min\n\n',...
        0.01*round(100*toc(stic)/60));
fclose(fid);


end




