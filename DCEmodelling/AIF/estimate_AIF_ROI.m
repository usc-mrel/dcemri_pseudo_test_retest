function [aifROI, vesselRegions, aifIdx] = estimate_AIF_ROI(img, opt)
%function aifROI = estimate_AIF_ROI(img, opt)
%   estimates the ROI for AIF estimation based on features (see below)
%
%
%   Inputs:
%
%       img             (complex) image data [Nx Ny Nz Nt]
%       opt             options
%         .size                         vector of dimension: [Nx Ny Nz Nt]
%         .SIdirection                  dimension of S/I direction
%         .ibolus                       time index of bolus arrival
%         .time                         time vector [Nt]
%         .voxel_size                   voxel size in mm
%         .mode2D                       (optional) if 3D data is not availble,
%                                       operate in 2D mode (not maintained)
%         .strip_transverse_vessels     (deprecated)
%         .erode_vessel                 (optional) apply vessel erosion
%         .threshold                    (optional) threshold for each vessel feature
%                                           .TTP                minimum time to peak in fraction of observation window
%                                           .FW80M              maximum full-width-80%-maximum
%                                           .peaktobaseline     minimum peak to baseline difference
%                                           .vesselness         minimum vesselness         
%
%   Outputs:
%       AIFROI          linear indices of AIF ROI
%       vesselRegions   cell strucutre of all vessels that have been
%                       detected
%       aifIdx          index of vessel in vesselRegions that is used as
%                       AIF
%
%   Example setup:
%       opt.fid                 = 1;                % log file handle
%       opt.size                = [256 256 160 50];
%       opt.SIdirection         = 1;
%       opt.ibolus              = 1;
%       opt.time                = 0:5:250;
%       opt.voxel_size          = [1 1 1];
%       aifROI = estimate_AIF_ROI(img, opt);
%
%
%   References:
%      [1] D. D. Shanbhag, S. N. Gupta, K. T. Rajamani, Y. Zhu, and R. Mullick, ?A generalized methodology for detection of vascular input function with dynamic contrast enhanced perfusion data,? in Proc. Intl. Soc. Mag. Reson. Med., vol. 20, no. 1, 2012, p. 3524.
%      [2] S. L. S. Chan and Y. Gal, ?Automatic Detection of Arterial Voxels in Dynamic Contrast-Enhanced MR Images of the Brain,? in International Conference on Digital Image Computing Techniques and Applications (DICTA), 2012, pp. 1?7.
%      [3] M. Jacobs, M. Benovoy, L.-C. Chang, A. E. Arai, and L.-Y. Hsu, ?Evaluation of an automated method for arterial input function detection for first-pass myocardial perfusion cardiovascular magnetic resonance,? Journal of Cardiovascular Magnetic Resonance, vol. 18, no. 1, p. 17, dec 2016.
%      [4] Frangi AF, Niessen WJ, Vincken KL, Viergever MA. Multiscale vessel enhancement filtering. 1998;1496:130?137 doi: 10.1007/BFb0056195.
%      [5] Chan SLS, Gal Y. Automatic Detection of Arterial Voxels in Dynamic Contrast-Enhanced MR Images of the Brain. In: International Conference on Digital Image Computing Techniques and Applications (DICTA). ; 2012. pp. 1?7. doi: 10.1109/DICTA.2012.6411710.
%      ... and plenty more
%
%   Yannick Bliesener, bliesene@usc.edu, 2018/2019
%
%

flag_erode_vessel                    = true;
flag_strip_transverse_vessels        = true;

%% put SI direction on third
perm = circshift(1:3, [0, 3-opt.SIdirection]);
perm = [perm 4];

img = permute(img, perm);

[Nx,Ny,Nz,Nt] = size(img);
N = Nx*Ny*Nz;

mode3D = (Nz > 3);

if isfield(opt, 'mode2D')
    mode3D = ~opt.mode2D;
end
if isfield(opt, 'strip_transverse_vessels')
    flag_strip_transverse_vessels = opt.strip_transverse_vessels;
end
if isfield(opt, 'erode_vessel')
    flag_erode_vessel = opt.erode_vessel;
end
if isfield(opt, 'voxel_size')
    voxel_size = opt.voxel_size(perm(1:3));
else
    voxel_size = [1 1 1];
end

%% set default threshold parameters

% minimum time to peak in fraction of observation window
threshold.TTP = 0.20;
% maximum full-width-80%-maximum
threshold.FW80M = 0.3;
% minimum peak to baseline difference
threshold.peaktobaseline = 0.2;
% minimum vesselness
threshold.vesselness = 0.001;

% threshold.minOverlap = 2/3;

% overwrite default parameters
if isfield(opt, 'threshold')
    for fn = fieldnames(opt.threshold)'
        threshold.(fn{1}) = opt.threshold.(fn{1});
    end
end

%% Compute relevant deatures for AIF detection

% normalize data
imgNormalized = img;
imgNormalized = imgNormalized ./ max(abs(imgNormalized(:)));

% baseline
baseline = mean(imgNormalized(:,:,:,1:opt.ibolus),4);

imgNormalized = imgNormalized - repmat(baseline, [1 1 1 opt.size(4)]);
imgNormalized = abs(imgNormalized);

% enhancement
% dImg = imgNormalized - mean(imgNormalized(:,:,:,1:opt.ibolus),4);

% peak value
[features.peakValue, features.peakTime] = max( imgNormalized, [], 4);

% subtraction masks
features.subMaskPeak2Last = features.peakValue - imgNormalized(:,:,:,end);
% features.subMaskPeak2First = features.peakValue - imgNormalized(:,:,:,1);
features.subMaskLast2First = mean(imgNormalized(:,:,:,end-3:end),4) - mean(imgNormalized(:,:,:,1:max(opt.ibolus-1,1)),4);
features.subMaskPeak2PreBolus = features.peakValue;
    % since the average of the pre-bolus is substracted this is zero

% vesselness
MRAimg = features.peakValue - imgNormalized(:,:,:,end);

% Note: How to interpret the sigmas?
% - It's a multiscale approach and each row corresponds to a scale
% - mutliply the sigmas with 6 (99% neighborhood) and you get roughly the size of the vessel
%   (in mm) that you want the filter to detect 
vesselSigmas = [
    0.5, 0.5, 0.5;
    1 1 1;
    2 2 2;
    5 5 5;
];
vesselSigmas = vesselSigmas ./ repmat(reshape(voxel_size, 1, 3), [size(vesselSigmas,1), 1]); % guess what: Matlab 2015 compatibility.....
features.vesselness = vesselness(MRAimg, vesselSigmas);

% time to peak
[~,I] = max( imgNormalized, [], 4);
features.TTP = (opt.time(I) - opt.time(opt.ibolus)) / (opt.time(Nt) - opt.time(1));
% features.TTP( features.TTP < 0 ) = 0;
clear I

% FW-80%-Maximum
z = imgNormalized > 0.8*repmat(features.peakValue, [1 1 1 Nt]);
z = reshape(z, N, Nt);
features.FW80M = -ones(Nx,Ny,Nz);
for i=1:N
    zrise = find(z(i,:),1,'first');
    if isempty(zrise) 
        zrise = 1;
    end
    zfall = find((~z(i,:))&((1:Nt)>zrise),1,'first');
    if isempty(zfall)
        zfall = Nt;
    end
    features.FW80M(i) = opt.time(zfall) - opt.time(zrise);
    features.FW80M(i) = features.FW80M(i) / (opt.time(Nt) - opt.time(1));
end
clear z

% z = imgNormalized > 0.9*repmat(features.peakValue, [1 1 1 Nt]);
% z = reshape(z, N, Nt);
% features.FW90M = -ones(Nx,Ny,Nz);
% for i=1:N
%     features.FW90M(i) = find(z(i,:),1,'last') - find(z(i,:),1,'first');
%     features.FW90M(i) = features.FW90M(i) / Nt;
% end
% clear z

% % AUC
% features.AUC = sum(dImg,4);
% 
% % normalize wrt AUC
% imgNormalizedAUC = dImg ./ repmat(features.AUC, [1 1 1 Nt]);
% 
% % first moment
% features.firstMoment = sum( imgNormalizedAUC .* repmat(reshape(opt.time, [1 1 1 Nt]), [Nx Ny Nz 1]), 4);

clear imgNormalized imgNormalizedAUC baseline

%% perform preselection for AIF detection

% strip background and non-enhancing pixel
Ipreselect = (features.subMaskPeak2Last > 0.05*max(features.subMaskPeak2Last(:)) );
% Ipreselect = (features.subMaskPeak2Last > 0.2*max(features.subMaskPeak2Last(:)) );

Iselect = Ipreselect & (features.TTP < threshold.TTP) & ...
    (features.FW80M < threshold.FW80M) & ...
    (features.subMaskPeak2PreBolus > threshold.peaktobaseline) & ...
    (features.vesselness > threshold.vesselness);

%% perform clean up and choose the vessel for AIF detection

if mode3D   % this is for 3D mode
    
    % strip isolated voxels
    for i = 1:Nz
        Iselect(:,:,i) = bwareaopen(Iselect(:,:,i),4);
    end
    clear i

    if sum( Iselect(:) ) < 1
        err_msg = 'Could not find AIF vessel';
        fprintf(opt.fid, '\t\t\t ERROR: %s \n', err_msg);
        error(err_msg)
    end

    % find volumes of connected components
    
    % this is for Matlab > R2017b
%     vessels = regionprops3(Iselect, 'VoxelIdxList', 'Volume', 'EigenVectors', 'EigenValues');
    
    % this is for Matlab < R2017b
    vessels = my_regionprops3(Iselect, [1 1 1], 'VoxelIdxList', 'Volume', 'EigenVectors', 'EigenValues', 'Zextend');

    % pick top volume vessels
%     [~, maxVolumeVessels] = maxk(vessels.Volume,top_n_vessel);
%     vessels = vessels(maxVolumeVessels,:);
%     [~, Izextend] = sort(vessels.Zextend,'descend');
%     
%     num_vessels = size(vessels,1);
%     
%     % retrieve principle axes
%     principle_axes = cell(1,num_vessels);
%     for nvessel=1:num_vessels
%         principle_axes{nvessel} = vessels.EigenVectors{nvessel}(:,1) * vessels.EigenValues{nvessel}(1);
%     end
%     principle_axes = cell2mat(principle_axes);
%     
%     % project onto positive SI direction
%     principle_axes = principle_axes(3,:);
%     
%     % pick vessel that is most aligned with B0
%     [~, maxParallelVessel] = max(principle_axes);

%     voxelList = vessels.VoxelIdxList{maxParallelVessel};
    
    [~, Isort] = sort(vessels.Volume, 'descend');
    vesselRegions = vessels.VoxelIdxList;
    vesselRegions = vesselRegions(Isort);
    
    % choose biggest vessel
    aifIdx = 1;
    voxelList = vesselRegions{aifIdx};
    
%     aifIdx = find(Isort == maxParallelVessel);

    Iselect = false(size(Iselect));
    Iselect(voxelList) = true;
    
    %% this section should be able to peel off the top part of the SS which is mostly transversal
    
%     if flag_strip_transverse_end_of_vessel
%         Ipreselect = Iselect;
% 
%         orient = zeros([Nz,1]);
%         for i = 1:Nz
% 
%             Ipreselect(:,:,i) = false;
% %             tmp_vessels = regionprops3(Ipreselect, 'Volume', 'EigenVectors');
%             tmp_vessels = my_regionprops3(Ipreselect, [1 1 1], 'Volume', 'EigenVectors');
% 
%             % pick top volume vessels
%             [~, tmp_maxVolumeVessels] = max([tmp_vessels.Volume]);
%             tmp_vessels = tmp_vessels(tmp_maxVolumeVessels,:);
%             num_vessels = size(tmp_vessels,1);
% 
%             if num_vessels < 1
%                 break
%             end
% 
%             % retrieve principle axes
%             tmp_principle_axes = cell(1,num_vessels);
%             for nvessel=1:num_vessels
%                 tmp_principle_axes{nvessel} = tmp_vessels.EigenVectors{nvessel}(:,1);
%             end
%             tmp_principle_axes = cell2mat(tmp_principle_axes);
% 
%             % project onto positive SI direction
%             orient(i) = tmp_principle_axes(3,:);
% 
%             if (i > 1) && (orient(i-1) > orient(i))
%                 break
%             end        
%         end
% 
% 
%         [maxOrientation,I] = max(orient);
%         I = find((orient) > 0.98*maxOrientation, 1, 'first');
%         Iselect(:,:,1:(I-1)) = 0;
% 
%     %     tmp_vessels = regionprops3(Iselect, 'VoxelIdxList', 'Volume');
%         tmp_vessels = my_regionprops3(Iselect, [1 1 1], 'VoxelIdxList', 'Volume');
% 
%         % pick top volume vessels
%         [~, tmp_maxVolumeVessels] = max([tmp_vessels.Volume]);
%         tmp_vessels = tmp_vessels(tmp_maxVolumeVessels,:);
%         num_vessels = size(tmp_vessels,1);
% 
%         Iselect = false(size(Iselect));
%         Iselect( tmp_vessels.VoxelIdxList{1} ) = true;
% 
%         clear tmp_* I
%     end
    %% Erode to exclude partial volume effect
    Ipreselect = Iselect;
    try
        SE = strel('sphere',1);
    catch e
        % Matlab R2015b compatibility
        load('vesselSE.mat', 'SE')
    end
    Iselect = imerode(Iselect, SE);
    
    if sum( Iselect(:) ) < 1
        err_msg = 'Could not find AIF vessel after PVE correction!';
        fprintf(opt.fid, '\t\t\t ERROR: %s \n', err_msg);
        error(err_msg)
    end
    
else    % 2D mode (mostly for debugging when only a single slice is chosen)
   
    Ipreselect = Iselect;
    
	% This thing is a cascade of erosions which are automatically tamed
    % down if found to agreesive to yield a good vessel!
    createMask = true;
    while createMask
        
        Iselect = Ipreselect;

        % erode bitmask to exclude PVE
        if flag_erode_vessel
            SE = strel('disk',1,4);
            for i = 1:Nz
                Iselect(:,:,i) = imerode(Iselect(:,:,i),SE);
            end
            clear SE i
        end

        % strip isolated voxels
        for i = 1:Nz
            Iselect(:,:,i) = bwareaopen(Iselect(:,:,i),4);
        end
        clear i 

        if flag_strip_transverse_vessels
            % stip lengthy objects
            % these are transverse vessels which do not fit the assumption of vessel || B0 
            % and hence do not fit the susceptibility model 
            for i = 1:Nz
                cc = bwconncomp(squeeze(Iselect(:,:,i))); 
                if cc.NumObjects > 0
                    stats = regionprops(cc, 'Eccentricity');

                    for k = 1:length(stats)
                       if stats(k).Eccentricity > 0.95

                           z = Iselect(:,:,i);
                           z( cc.PixelIdxList{k}(:) ) = 0;
                           Iselect(:,:,i) = z;

                       end
                    end
                end
            end
            clear cc z i k stats
        end
        
        if sum( Iselect(:) ) < 1
            if flag_erode_vessel
                % retry without eroding
                flag_erode_vessel = false;
                continue
            else
                if flag_strip_transverse_vessels
                    % retry whith keeping transverse vessels
                    % OBS: this will violate the modelling assumption in
                    % the AIF estimation suscpetibity kernel but its
                    % probably better than failure!
                    flag_strip_transverse_vessels = false;
                    continue
                else
                    createMask = false;
                    break
                end
            end
        else
            createMask = false;
            break
        end
        
    end
        
    if sum( Iselect(:) ) < 1
        err_msg = 'Could not find AIF vessel';
        fprintf(opt.fid, '\t\t\t ERROR: %s \n', err_msg);
        error(err_msg)
    end
    
	% find areas of connected components
    vesselArea          = cell(1,Nz);
    vesselRegions       = cell(1,Nz);
    Nvessels = 0;
    for slice = 1:Nz
        stats = regionprops(Iselect(:,:,slice), 'PixelList', 'Area');
        vesselArea{slice}          = { stats.Area };
        vesselRegions{slice}       = { stats.PixelList };
        clear stats
    
        Nvessels = Nvessels+length(vesselRegions{slice});
    end
    clear slice
    
    
    % convert to 3D coordinates
    for slice = 1:Nz
        for ivessel = 1:length(vesselRegions{slice})
            curPixelList = [vesselRegions{slice}{ivessel}(:,2), vesselRegions{slice}{ivessel}(:,1), slice*ones(size(vesselRegions{slice}{ivessel}(:,1)))];
            vesselRegions{slice}{ivessel} = sub2ind([Nx Ny Nz], curPixelList(:,1), curPixelList(:,2), curPixelList(:,3));
        end
    end
    
    vesselArea    = horzcat(vesselArea{:});
    vesselRegions = horzcat(vesselRegions{:}); 

    % find vessel with maximum area
    [~, Isort] = sort(cell2mat(vesselArea), 'descend');
    vesselArea    = vesselArea(Isort);
    vesselRegions = vesselRegions(Isort);
    
    voxelList = vesselRegions{1};
    
	Iselect = false(size(Iselect));
    Iselect(voxelList) = true;
    
end


% figure(666),
% imgt = abs(img(:,:,:,opt.ibolus+2));
% red = cat(4, ones(Nx,Ny,Nz), zeros(Nx,Ny,Nz), zeros(Nx,Ny,Nz));
% ROI = Iselect;
% 
% imgt = reshape(imgt, [Nx Ny 2 3]);
% imgt = permute(imgt, [1 3 2 4]);
% imgt = reshape(imgt, [Nx*2 Ny*3]);
% 
% red = reshape(red, [Nx Ny 2 3 3]);
% red = permute(red, [1 3 2 4 5 ]);
% red = reshape(red, [Nx*2 Ny*3 3]);
% 
% ROI = reshape(ROI, [Nx Ny 2 3]);
% ROI = permute(ROI, [1 3 2 4]);
% ROI = reshape(ROI, [Nx*2 Ny*3]);
% 
% imshow(imgt, [])
% axis image off, hold on
% h = imshow( red );
% set(h, 'AlphaData', 0.6*ROI)

% figure(666),
% for i=1:Nz
%     subplot(6,Nz/6,i), imshow(abs(img(:,:,i,opt.ibolus+2)), [])
%     axis equal off, hold on
%     red = cat(3, ones(Nx,Ny), zeros(Nx,Ny), zeros(Nx,Ny));
% 
%     h = imshow( red );
%     set(h, 'AlphaData', 0.6*Iselect(:,:,i))
% end


Iselect = ipermute(Iselect, perm(1:3));
aifROI = find(Iselect);

for v = 1:length(vesselRegions)
    vessel = false(Nx,Ny,Nz);
    vessel(vesselRegions{v}) = true;
    vessel = ipermute(vessel, perm(1:3));
    vesselRegions{v} = find(vessel(:));
end

end