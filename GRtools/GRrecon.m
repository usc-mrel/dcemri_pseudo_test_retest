function [imginf, imgPRE, imgPEAK, imgPOST, imgVS, imgR, OUT] = GRrecon(GRopt,path, outputfile_suffix)
% function [imginf, imgPRE, imgPEAK, imgPOST, imgVS, imgR, OUT] = GRrecon(GRopt,path, outputfile_suffix)
%
%   Function to reconstruct dynamic undersampled data
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
%       imginf              all data merged into one time frame
%       imgPRE              pre contrast anatomic image
%       imgPEAK             peak contrast anatomic image
%       imgPOST             post contrast anatomic image
%       imgVS               view-shared
%       imgR                anatomic image time series
%       OUT                 additional output information
%
%   Usage:
%       GRopt = GRoptset(64, 'perm');
%       GRopt.mtrxH = [256 240 120];
%       GRopt.DO_SPSENSE = 1;
%       [~, ~, ~, ~, imgR, OUT] = GRrecon(GRopt,[]);
%
%
%   Author
%       R. Marc Lebel, mlebel@gmail.com
%       Yannick Bliesener, bliesene@usc.edu, MREL at USC, 2019/2020
%       
%

%   Start clock and open log file
stic = tic;
fid = fopen(['GRlog_' GRopt.app '.txt'],'w');
fprintf(fid,'--- Cartesian radial image reconstruction ---\n');
fprintf(fid,'R. Marc Lebel\nmlebel@gmail.com\n\n');

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

%   Remove dummy pulses
% NOTE: STARDCE probably won't have that
if hdr.ImageData.user46
    fprintf(fid,'Removing dummy pulses...\n');tic;
    NSTST = hdr.ImageData.user46;
    k = k(:,(NSTST+1):end,:,:,:);
    petab = petab((NSTST+1):end,:);
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
    NCAL = 0;
end

%   Strip T1 mapping
if hdr.ImageData.user30    
    nFA = hdr.ImageData.user29;
    nTR = nFA * floor(hdr.ImageData.user30/nFA);
    
    % just overwrite the defaul values
    GRopt.VFA.nFA = nFA;
    GRopt.VFA.nTR = nTR;
else
    warning('Assuming default VFA information!')
    fprintf(fid,'Could not find VFA information!\n');
    fprintf(fid,'->Using default VFA values!\n');
    
    nFA = GRopt.VFA.nFA;
    % number of flip angle prior to DCE acquisition flip angle
    nTR = GRopt.VFA.nTR;
    % number of TRs for VFA
end
% stripping VFA, coil sensitivty calibration region, eddy current
% calibration region, and dummy pulses
indstart = find(tstamp>(nTR + NCAL + NEDDY + NSTST)*tr,1);
fprintf(fid,'Removing T1 mapping data (0 - %3.2f sec)\n',tstamp(indstart));
k = k(:,indstart:end,:,:,:);
petab = petab(indstart:end,:);
tstamp = tstamp(indstart:end);
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

% create output file
if nargin < 3 || isempty(outputfile_suffix)
    matname = ['GRrecon_' GRopt.app '.mat'];
else
    matname = ['GRrecon_' GRopt.app '_' outputfile_suffix '.mat'];
end
save(matname,'OUT', 'githash', matfile_version);


%   Estimate coil sensitvities
fprintf(fid,'Estimating coil sensitivities...');tic
[S,sf] = GRcoilsens(kS,NCAL);
k = sf*k;
clear kS
fprintf(fid,'done (%g sec)\n',toc);


%   Check for scaling factor from prior T1/Mo mapping
sf2 = 8192;
if exist('GRrecon_T1Mo.mat','file')
    load('GRrecon_T1Mo.mat','sf2');
end
k = sf2 * k;
sf2 = 1;    % important to prevent it from being reapplied during phase correction

% for debugging only: downsample
if GRopt.debug > 0
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
    pfle.xRes = length(GRopt.debug);
    np = length(GRopt.debug);
    GRopt.mtrxH(1) = length(GRopt.debug);
    clear ind tmp_opt
    ind = GRopt.debug;
end

%   Open parallel pool
ncores = feature('numCores');
parpl = gcp('nocreate');
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


%   Get generic SparseSENSE recon options
opt = SPSENSE_optset;
opt.FDorder1 = {1,2,3};
opt.Dobj     = GRopt.Dobj;
opt.MaxIter  = GRopt.SP.MaxIter;
opt.MinIter  = 6;
opt.class    = GRopt.class;
opt.gpuid    = GRopt.gpuid;


%% %%%% STATIC IMAGE RECON %%%%%%
imginf = [];
if GRopt.DO_STATIC
    fprintf(fid,'\n\n---- STATIC IMAGE RECON (ITERATIVE) ----\n');
    imginf = STrecon(k,petab,pfle,tstamp);
    
    %   Save
    fprintf(fid,'Phase correcting...');tic;
    imgInf = PHcorrection(imginf,sf2);
    %imgInf = flip(imgInf,3);    %   Empirical
    fprintf(fid,'done (%g sec)\n',toc);
    fprintf(fid,'Saving .mat...');tic;
    save(matname,'imgInf','-append');
    fprintf(fid,'done (%g sec)\n',toc);
    if GRopt.DCM
        fprintf(fid,'Writing DICOMs...');tic;
        GRwrite_dicom(imgInf,pfle,hdr,'ImgINF','ImgINF','Entire scan',0,GRopt.DCM_install);
        fprintf(fid,'done (%g sec)\n',toc);
    end
end
%%%%%% STATIC IMAGE RECON %%%%%%




%% %%%% PRE/PEAK/POST IMAGE RECON %%%%%%
imgPRE  = [];
imgPEAK = [];
imgPOST = [];
if GRopt.DO_PRE_POST
    fprintf(fid,'\n\n---- PRE/PEAK/POST CONTRAST IMAGE RECON (ITERATIVE) ----\n');
    opt.IREF = imginf;
    opt.img0 = imginf;
    [imgPRE,imgPEAK,imgPOST] = PPrecon(k,petab,pfle,tstamp);
    
    %   Save
    fprintf(fid,'Phase correcting...');tic;
    imgPRE  = PHcorrection(imgPRE,sf2);
    imgPEAK = PHcorrection(imgPEAK,sf2);
    imgPOST = PHcorrection(imgPOST,sf2);
    %imgPRE = flip(imgPRE,3);
    %imgPOST = flip(imgPOST,3);
    fprintf(fid,'done (%g sec)\n',toc);
    fprintf(fid,'Saving .mat...');tic;    
    save(matname,'imgPRE','imgPEAK','imgPOST','-append');
    fprintf(fid,'done (%g sec)\n',toc);
    if GRopt.DCM
        fprintf(fid,'Writing DICOMs...');tic;
        GRwrite_dicom(imgPRE,pfle,hdr,'ImgPRE','imgPRE','Pre-contrast',1,GRopt.DCM_install);
        GRwrite_dicom(imgPEAK,pfle,hdr,'ImgPEAK','imgPEAK','Peak-contrast',2,GRopt.DCM_install);
        GRwrite_dicom(imgPOST,pfle,hdr,'ImgPOST','imgPOST','Post-contrast',3,GRopt.DCM_install);
        fprintf(fid,'done (%g sec)\n',toc);
    end
end
%%%%%% PRE/PEAK/POST IMAGE RECON %%%%%%

%% %%%% VIEW-SHARING RECON %%%%%%
imgVS = [];
if GRopt.DO_VIEW_SHARING
    fprintf(fid,'\n\n---- VIEW-SHARING RECON ----\n');
    imgvs = VSrecon(k,imginf,petab,pfle,OUT.binE,tstamp);
    
    %   Save
    fprintf(fid,'Phase correcting...');tic;
    imgVS = PHcorrection(imgvs,sf2);
    %imgVS = flip(imgVS,3);
    fprintf(fid,'done (%g sec)\n',toc);
    fprintf(fid,'Saving .mat...');tic;
    save(matname,'imgVS','-append');
    fprintf(fid,'done (%g sec)\n',toc);
    if GRopt.DCM && false
        fprintf(fid,'Writing DICOMs...');tic;
        GRwrite_dicom(imgVS,pfle,hdr,'ImgVS','imgVS','View-sharing',3,GRopt.DCM_install);
        fprintf(fid,'done (%g sec)\n',toc);
    end
end
%%%%%% VIEW-SHARING RECON %%%%%%




%% %%%% DYNAMIC IMAGE RECON %%%%%%
imgR = [];
if GRopt.DO_SPSENSE
    fprintf(fid,'\n\n---- DYNAMIC IMAGE RECON (ITERATIVE) ----\n');
    opt.S = S;
    
    %   EXPERIMENTAL
    if exist('imgvs','var') && ~isempty(imgvs)
        opt.IREF = svd_smooth(imgvs,opt.Npc);
        opt.img0 = opt.IREF;
    else
        opt.IREF = [];
        opt.img0 = [];
    end
    clear imgvs
    
    %   Recon image
    imgR = DYrecon(k,petab,pfle,OUT.binE,tstamp);
    
    %   Save
    try
        fprintf(fid,'Saving .mat...');tic;
        save(matname,'imgR','-append');
        fprintf(fid,'done (%g sec)\n',toc);
        fprintf(fid,'Phase correcting...');tic;
        imgR = PHcorrection(imgR,sf2);
        %imgR = flip(imgR,3);
        fprintf(fid,'done (%g sec)\n',toc);
        if GRopt.DCM
            fprintf(fid,'Writing DICOMs...');tic;
            GRwrite_dicom(imgR,pfle,hdr,'ImgDY','imgDY','SparseSENSE dynamic',4,GRopt.DCM_install);
            fprintf(fid,'done (%g sec)\n',toc);
        end
    catch e
       save('backup.mat', '-v7.3')
    end
end
%%%%%% DYNAMIC IMAGE RECON %%%%%%




%% %%%% ANGIO %%%%%%
if GRopt.DO_ANGIO
    fprintf(fid,'\n\n---- ANGIOGRAPHY ----\n');
    fprintf(fid,'Loading images...');tic;
    load('GRrecon_angio','imgPRE','imgR');
    fprintf(fid,'done (%g sec)\n',toc);
    
    fprintf(fid,'Computing projections...');tic
    make_angio(imgPRE,imgR,OUT,GRopt.DCM);s
    fprintf(fid,'done (%g sec)\n',toc);
end
%%%%%% ANGIO %%%%%%




%   End
save(matname,'OUT','GRopt','-append');
fprintf(fid,'\n ---- DONE ----\nReconstruction time: %g min\n\n',...
        0.01*round(100*toc(stic)/60));
fclose(fid);


if ~isempty(parpl)
    %   Close parallel pool
    delete(parpl);
    pc.Jobs.delete;
end


%% %%%% NESTED FUNCTIONS %%%%%%

    function stimg = STrecon(k_in,pe_in,pfle_in,tstamp_in) 

        %   Recon single time frame data
        fprintf(fid,'Reshaping raw data...');tic;
        [kst,Ust,kwst] = GRreshape2(k_in,pfle_in,pe_in,Inf,tstamp_in,0);
        kst = iFastFT(kst,1,1);
        kst = fftshift(kst,2);
        kst = fftshift(kst,3);
        Ust = fftshift(fftshift(Ust,2),3);
        kwst = fftshift(fftshift(kwst,2),3);
        fprintf(fid,'done (%g sec)\n',toc);
        
         %%%   Set penatlty   %%%
        opt.lambda = GRopt.ST.lambda;
        opt.compression_levels = GRopt.ST.compression_levels;
        OUT.lambdaST = opt.lambda;
        OUT.compression_levelsST = opt.compression_levels;
        fprintf(fid,'Lambda = [');
        fprintf(fid,' %g',opt.lambda);
        fprintf(fid,']\n');
        
        
        %   Recon full size, single time frame image
        fprintf(fid,'Reconstructing full size image...');tic;
        opt.U = find(Ust);
        R = nv*ns/sum(sum(squeeze(Ust(1,:,:,1,1)),1),2);clear Ust;
        kpack = kst(opt.U);clear kst;
        kpackwt = kwst(opt.U);clear kwst;
%         kpackwt = ones(size(kpackwt),class(kpack));    %%% TODO %%%
%         opt.kW = kpackwt(:);
        opt.kW = ones(1,GRopt.class);
        opt.size = [np nv ns 1 nr];
        opt.S = S;
        [stimg,opttmp] = SPSENSE_recon(kpack,opt);
        fprintf(fid,'done (%g sec; %d iter)\n',toc,opttmp.RI.Iter);
        fprintf(fid,'Acceleration:      %gx\n',R);
        if opt.lambda(2) > 0
            fprintf(fid,'Wavelet transform: ');
            fprintf(fid,'{[');
            for iind = 1:length(opttmp.worder)
                fprintf(fid,' %d',opttmp.worder{iind});
                if iind < length(opttmp.worder)
                    fprintf(fid,'],[');
                end
            end
            fprintf(fid,']}\n');
        end
        clear opttmp kpack iind kpackwt;
        %%%%%% DONE STATIC IMAGE %%%%%%
    end


    function [preimg, peakimg, postimg] = PPrecon(k_in,pe_in,pfle_in,tstamp_in)

        %   Recon single time frame data
        fprintf(fid,'Reshaping raw data...');tic;
        bins = [GRopt.PP.windows(1) 0.5*(tstamp_in(end)+GRopt.PP.windows(2)) tstamp_in(end)];
            % this should be
            %   bolus arrival
            %   mid point between bolus arrival and end point
            %   end point
        [kpre,Upre,kwpre] = GRreshape2(k_in,pfle_in,pe_in,bins,tstamp_in,0);
        [kpeak,Upeak,kwpeak] = GRreshape2(k_in,pfle_in,pe_in,bins,tstamp_in,0);
        [kpst,Upst,kwpst] = GRreshape2(k_in,pfle_in,pe_in,bins,tstamp_in,0);
        kpre = kpre(:,:,:,1,:);
        kpeak = kpeak(:,:,:,2,:);
        kpst = kpst(:,:,:,3,:);
        Upre = Upre(:,:,:,1,:);
        Upeak = Upeak(:,:,:,2,:);
        Upst = Upst(:,:,:,3,:);
        kwpre= kwpre(:,:,:,1,:);
        kwpeak= kwpeak(:,:,:,2,:);
        kwpst= kwpst(:,:,:,3,:);
        kpre = iFastFT(kpre,1,1);
        kpre = fftshift(kpre,2);
        kpre = fftshift(kpre,3);
        kpeak = iFastFT(kpeak,1,1);
        kpeak= fftshift(kpeak,2);
        kpeak = fftshift(kpeak,3);
        kpst = iFastFT(kpst,1,1);
        kpst = fftshift(kpst,2);
        kpst = fftshift(kpst,3);
        Upre = fftshift(fftshift(Upre,2),3);
        Upeak = fftshift(fftshift(Upeak,2),3);
        Upst = fftshift(fftshift(Upst,2),3);
        kwpre = fftshift(fftshift(kwpre,2),3);
        kwpeak = fftshift(fftshift(kwpeak,2),3);
        kwpst = fftshift(fftshift(kwpst,2),3);
        fprintf(fid,'done (%g sec)\n',toc);
        
        %   Set options
        opt.S = S;
        opt.size = [np nv ns 1 nr];
        opt.lambda = GRopt.PP.lambda;
        opt.compression_levels = GRopt.PP.compression_levels;
        fprintf(fid,'Lambda = [');
        fprintf(fid,' %g',opt.lambda);
        fprintf(fid,']\n');
        OUT.lambdaPP = opt.lambda;
        OUT.compression_levelsPP = opt.compression_levels;
        
        
        %   Recon full size pre-contrast image
        fprintf(fid,'Reconstructing Pre-contrast image...');tic;
        opt.U = find(Upre);
        R = nv*ns/sum(sum(squeeze(Upre(1,:,:,1,1)),1),2);
        kpack = kpre(opt.U);clear kpre Upre;
%         kpackwt = kwpre(opt.U);
%         opt.kW = kpackwt(:);
        clear kwpre;
        opt.kW = 1;
        [preimg,opttmp] = SPSENSE_recon(kpack,opt);
        fprintf(fid,'done (%g sec; %d iter)\n',toc,opttmp.RI.Iter);
        fprintf(fid,'Acceleration:      %gx\n',R);
        
        %   Recon full size peak-contrast image
        fprintf(fid,'Reconstructing Peak-contrast image...');tic;
        opt.U = find(Upeak);
        R = nv*ns/sum(sum(squeeze(Upeak(1,:,:,1,1)),1),2);
        kpack = kpeak(opt.U);clear kpeak Upeak;
        %         kpackwt = kwpre(opt.U);
        %         opt.kW = kpackwt(:);
        clear kwpeak;
        opt.kW = 1;
        [peakimg,opttmp] = SPSENSE_recon(kpack,opt);
        fprintf(fid,'done (%g sec; %d iter)\n',toc,opttmp.RI.Iter);
        fprintf(fid,'Acceleration:      %gx\n',R);
        
        %   Recon full size post-contrast image
        fprintf(fid,'Reconstructing Post-contrast image...');tic;
        opt.U = find(Upst);
        R = nv*ns/sum(sum(squeeze(Upst(1,:,:,1,1)),1),2);
        kpack = kpst(opt.U);clear kpst Upst;
%         kpackwt = kwpst(opt.U);
%         opt.kW = kpackwt(:);
        clear kwpst;
        [postimg,opttmp] = SPSENSE_recon(kpack,opt);
        fprintf(fid,'done (%g sec; %d iter)\n',toc,opttmp.RI.Iter);
        fprintf(fid,'Acceleration:      %gx\n',R);
        clear opttmp kpack iind kpackwt;
        %%%%%% DONE PRE/POST CONTRAST IMAGE %%%%%%
    end


    function vsimg = VSrecon(kvs,kfl,pe_in,pfle_in,binE_in,tstamp_in)
        
        %   Get key parameters
        [np_vs, nv_vs, ns_vs, nr_vs] = get_key_params2(pfle_in);
        
        fprintf(fid,'Reshaping raw data...');tic;
        kfl = fFastFT(fSENSE(kfl,opt),[1 2 3],1);
        [kvs,Uvs] = GRreshape2(kvs,pfle_in,pe_in,binE_in,tstamp_in,0);
        ph = size(kvs,4);
        fprintf(fid,'done (%g sec)\n',toc);
        
        fprintf(fid,'Convolving data...');tic;
        Uvs = find(Uvs(1,:,:,:,1));
        [yi,zi,ti] = ind2sub([nv_vs ns_vs ph],Uvs);
        clear Uvs;
        pet = [yi zi ti];
        window = ceil(15/median(diff(binE_in)));    %   +- 15 second window
        fwhm = GRopt.VS.fwhm;     %    1.25 frame FWHM for Gaussian weighting
        for iind = 1:nr_vs
            [kvs(:,:,:,:,iind),Uvs(:,:,:,:,iind)] = convVS(kvs(:,:,:,:,iind),pet,[np_vs nv_vs ns_vs ph 1],window,fwhm);
            fprintf(fid,' (%d)',iind);
        end
        clear pet yi zi ti fwhm window
        for iind = 1:ph
            tmp = kvs(:,:,:,iind,:);
            tmp2 = tmp == 0;
            tmp(tmp2) = kfl(tmp2);
            kvs(:,:,:,iind,:) = tmp;
        end
        clear kfl tmp tmp2
        fprintf(fid,' ...done (%g sec)\n',toc);
        
        fprintf(fid,'Converting to image space...');tic;
        for iind = 1:nr_vs
            kvs(:,:,:,:,iind) = iFastFT(kvs(:,:,:,:,iind),[1 2 3],1);
        end
        vsimg = kvs;
        clear kvs
        fprintf(fid,'done (%g sec)\n',toc);
        
        %   Combine coils
        fprintf(fid,'Combining array coil images...');tic;
        opt.size = size(vsimg);
        opt.S = S;
        vsimg = iSENSE(vsimg,opt);
        fprintf(fid,'done (%g sec)\n',toc);
        %%%%%% DONE VIEW-SHARING RECON %%%%%%
        
    end


    function imgDY = DYrecon(k,pe_in,pfle_in,binE_in,tstamp_in)
        
        %   Get key parameters
        [np_dy, nv_dy, ns_dy, nr_dy] = get_key_params2(pfle_in);
        
        %   Rebin data and shift
        fprintf(fid,'Reshaping raw data...');tic;
%         [kdy,Udy,kwdy] = GRreshape(k,pe_in,par_in,binE_in,tstamp_in,0);   %   TODO
        [kdy,Udy] = GRreshape2(k,pfle_in,pe_in,binE_in,tstamp_in,0);
        ph = size(kdy,4);
        R = nv_dy*ns_dy/median(sum(sum(squeeze(Udy(1,:,:,:,1)),1),2));
        for iind = 1:nr_dy
            kdy(:,:,:,:,iind)  = iFastFT(kdy(:,:,:,:,iind),1,1);
            kdy(:,:,:,:,iind)  = fftshift(fftshift(kdy(:,:,:,:,iind),2),3);
            Udy(:,:,:,:,iind)  = fftshift(fftshift(Udy(:,:,:,:,iind),2),3);
%             kwdy(:,:,:,:,iind) = fftshift(fftshift(kwdy(:,:,:,:,iind),2),3);
        end
        clear iind jind
        kdy = kdy(Udy);
        kdy = reshape(kdy,[np_dy numel(kdy)/np_dy]);
        
        %   Transer regularization factors
        opt.lambda = GRopt.SP.lambda;           %   Current best
        opt.compression_levels = GRopt.SP.compression_levels;
        OUT.lambdaDY = opt.lambda;
        OUT.compression_levelsDY = opt.compression_levels;
        
        %   Break into readout segments for recon
        if GRopt.SP.nbl > np_dy
            nbl = np_dy;
        else
            nbl = GRopt.SP.nbl;
        end
        
        %   Remove FD in RO direction if size of 1
        if np_dy/nbl == 1
            opt.FDorder1 = {2,3};
        end
        
        
        %   Reshape data for easier access to blocks
        kdy = reshape(kdy,[np_dy/nbl nbl numel(kdy)/np_dy]);
        kdy = permute(kdy,[1 3 2]);
        Udy = reshape(Udy,[np_dy/nbl nbl nv_dy ns_dy ph nr_dy]);
        Udy = permute(Udy,[1 3 4 5 6 2]);
        Udy = Udy(:,:,:,:,:,1); % assumes all blocks are the same
        Sdy = reshape(S,[np_dy/nbl nbl nv_dy ns_dy 1 nr_dy]);
        Sdy = permute(Sdy,[1 3 4 5 6 2]);
        if ~isempty(opt.img0)
            useI0 = true;   %   This variable avoids making I0 a braodcast variable in the parfor loop
            I0 = reshape(opt.img0,[np_dy/nbl nbl nv ns ph]);
            I0 = permute(I0,[1 3 4 5 2]);
        else
            useI0 = false;
            I0 = zeros([1 1 1 1 nbl],opt.class);    %   This allows the parfor to slice the variable
        end
        if ~isempty(opt.IREF)
            useIREF = true; %   This variable avoids making I0 a braodcast variable in the parfor loop
            IREF = reshape(opt.IREF,[np_dy/nbl nbl nv ns ph]);
            IREF = permute(IREF,[1 3 4 5 2]);
        else
            useIREF = false;
            IREF = zeros([1 1 1 1 nbl],opt.class);  %   This allows the parfor to slice the variable
        end
        
        
        %   Set some options and clean up others to remove large unused variables
        opt.size = [np_dy/nbl nv_dy ns_dy ph nr_dy];
        opt.U = find(Udy(:,:,:,:,:));
        opt.S = [];
        opt.IREF = [];
        opt.img0 = [];
        opt.kW = 1;
        fprintf(fid,'done (%g sec)\n',toc);
        fprintf(fid,'Lambda = [');
        fprintf(fid,' %g',opt.lambda);
        fprintf(fid,']\n');

        
        
        %   Recon subset of lines for PCA
        %   Not needed with block PCA - could remove
        stic1 = tic;
        if false && GRopt.SP.nbl>1 && GRopt.SP.lambda(6) && isempty(GRopt.SP.BF)
            fprintf(fid,'Estimating Principle Components...');
            bl_ind = round(linspace(1,np_dy,max(GRopt.parallel/2,1)+2));
            bl_ind = unique(bl_ind(2:end-1));
            kdy_sm = kdy(:,:,bl_ind);
            Sdy_sm = Sdy(:,:,:,:,:,bl_ind);
            imgDY = zeros([np_dy/nbl nv_dy ns_dy ph length(bl_ind)],opt.class);
            parfor (ibl = 1:length(bl_ind),max(GRopt.parallel/2,1))
                %for ibl = 1:length(bl_ind);
                kpack = kdy_sm(:,:,ibl);
                opt2 = opt;
                opt2.S = Sdy_sm(:,:,:,:,:,ibl);
                
                %   Recon readout block
                imgDY(:,:,:,:,ibl) = SPSENSE_recon(kpack,opt2);
                kpack = [];
            end
            imgDY = double(permute(imgDY,[1 5 2 3 4]));
            imgDY = reshape(imgDY,[length(bl_ind) nv ns ph]);
            [~,opt.BF] = svd_smooth(imgDY,opt.Npc,0);
            OUT.BF = opt.BF;
            GRopt.SP.BF = opt.BF;
            fprintf(fid,'done (%g sec)\n',toc);
        else
            OUT.BF = GRopt.SP.BF;
            opt.BF = GRopt.SP.BF;
        end
        
        if ~isempty(parpl)
            GRopt.parallel = min(nbl, parpl.NumWorkers);
        end        
        
        %   Loop through all readout segments
        fprintf(fid,'Reconstructing entire image...');
        bl_ind = 1:nbl;
        imgDY = zeros([np_dy/nbl nv_dy ns_dy ph nbl],opt.class);
        parfor (ibl = bl_ind,max(GRopt.parallel,1))
%         for ibl = bl_ind
            tic;
            kpack = kdy(:,:,ibl);
            opt2 = opt;
            opt2.S = Sdy(:,:,:,:,:,ibl);
            
            %   Recon readout block
            [imgDY(:,:,:,:,ibl),opttmp] = SPSENSE_recon(kpack,opt2);
            
            %   Save some info
            blTm(ibl) = toc;
            blIt(ibl) = opttmp.RI.Iter;
            kpack = []; opttmp = [];
        end
        for ibl = bl_ind
            fprintf(fid,'\n\tblock %d (%g sec; %d iter)',ibl,blTm(ibl),blIt(ibl));
        end
        clear npind nbl ibl opttmp kpack kpackw ibl opt2 blTm blIt
        
        
        %   Reshape image
        imgDY = permute(imgDY,[1 5 2 3 4]);
        imgDY = reshape(imgDY,[np nv ns ph]);
        
        
        fprintf(fid,'\n\tdone (%g sec)\n',toc(stic1));
        fprintf(fid,'Median acceleration: %gx\n',R);
        
    end


    function [im,ph] = PHcorrection(im,scle,ph)
        
        %   Create phase map, if not provided
        if nargin < 3 || isempty(ph)
            fwidth = min(round(0.8*[np nv ns]),160);
            H = FermiFilt([np nv ns],fwidth,1/20);
            
            ph = zeros(size(im),'like',im);
            parfor iind = 1:size(im,4)
                im2 = im(:,:,:,iind);
                im2 = fFastFT(im2,[1 2 3],1);
                im2 = im2 .* H;
                im2 = iFastFT(im2,[1 2 3],1);
                ph(:,:,:,iind) = exp(-1i*angle(im2));
            end
            clear H im2 fwidth
        end
        
        %   Make sure ph is the same size as the image
        if size(ph,4) == 1
            ph = repmat(ph,[1 1 1 size(im,4)]);
        end
        
        %   Apply phase correction
        %   Loop through time points
        parfor iind = 1:size(im,4)
            im(:,:,:,iind) = im(:,:,:,iind) .* ph(:,:,:,iind);
        end
        
        %   Take real part
        im = real(im);
        
        %   Scale
        im = scle*im;
        
    end

end


%%%%%% ACTUAL SUB FUNCTIONS %%%%%%
function H = FermiFilt(sz,N,T)

%   Dim 1
E = N(1)/sz(1);
x = -1:2/(sz(1)-1):1;
H1 = 1./(1+exp((abs(x)-E)/T))';
H = repmat(H1,[1,sz(2),sz(3)]);

%   Dim 2
E = N(2)/sz(2);
x = -1:2/(sz(2)-1):1;
H1 = reshape(1./(1+exp((abs(x)-E)/T)),[1 sz(2) 1]);
H = H.*repmat(H1,[sz(1),1,sz(3)]);

%   Dim3
E = N(3)/sz(3);
x = -1:2/(sz(3)-1):1;
H1 = reshape(1./(1+exp((abs(x)-E)/T)),[1 1 sz(3)]);
H = H.*repmat(H1,[sz(1),sz(2),1]);

end

function H = GaussFilt(sz,N)

%   Dim 1
sd = N(1)/sz(1);
x = -1:2/(sz(1)-1):1;
H1 = reshape(exp(-x.^2/(2*sd^2)),[sz(1) 1 1]);
H = repmat(H1,[1,sz(2),sz(3)]);

%   Dim 2
sd = N(2)/sz(2);
x = -1:2/(sz(2)-1):1;
H1 = reshape(exp(-x.^2/(2*sd^2)),[1 sz(2) 1]);
H = H.*repmat(H1,[sz(1),1,sz(3)]);

%   Dim3
sd = N(3)/sz(3);
x = -1:2/(sz(3)-1):1;
H1 = reshape(exp(-x.^2/(2*sd^2)),[1 1 sz(3)]);
H = H.*repmat(H1,[sz(1),sz(2),1]);

if ~isa(sz,'single')
    H = single(H);
end

end

%#ok<*NASGU>
%#ok<*UNRCH>
%#ok<*AGROW>
