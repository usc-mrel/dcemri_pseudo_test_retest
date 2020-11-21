function theta = estimateVesselTheta(vesselROI, aifROI, opt)
%function theta = estimateVesselTheta(vesselROI, opt)
%   Estimates the angle between the vessel and B0 based on the vessels
%   principle component
%
%
%
% Yannick 2020
%

if numel(vesselROI) == prod(opt.size(1:3))
    vesselROI = find(vesselROI);
end

vesselMask = false(opt.size(1:3));
vesselMask(vesselROI(:)) = true;

if isempty(aifROI)
   aifROI = find(vesselROI);
end

if numel(aifROI) == prod(opt.size(1:3))
    aifROI = find(aifROI);
end

Nx = opt.size(1);
Ny = opt.size(2);
Nz = opt.size(3);
voxel_size = opt.voxel_size;
SIdirection = opt.SIdirection;

theta = -inf(length(aifROI),1);

% this corresponds to 20mm isotropic voxel size centered at the respective
% voxel of interest
% neighborhoodSize
neighborhoodSize = ceil(20./opt.voxel_size / 2);

parfor i = 1:length(aifROI)
   
    [x, y, z] = ind2sub([Nx Ny Nz], aifROI(i));
    
    neighborhood = vesselMask(max(x-neighborhoodSize(1),1):min(x+neighborhoodSize(1),Nx),max(y-neighborhoodSize(2),1):min(y+neighborhoodSize(2),Ny), max(z-neighborhoodSize(3),1):min(z+neighborhoodSize(3),Nz));
    
    tmp_vessels = my_regionprops3(neighborhood, voxel_size, 'Volume', 'EigenVectors');
    
    % pick top volume vessels
    [~, tmp_maxVolumeVessels] = max([tmp_vessels.Volume]);
    tmp_vessels = tmp_vessels(tmp_maxVolumeVessels,:);
    num_vessels = size(tmp_vessels,1);
    
    if num_vessels < 1
        theta(i) = nan;
        continue
    end
    
    % retrieve principle axes
    tmp_principle_axes = cell(1,num_vessels);
    for nvessel=1:num_vessels
        tmp_principle_axes{nvessel} = tmp_vessels.EigenVectors{nvessel}(:,1);
    end
    tmp_principle_axes = cell2mat(tmp_principle_axes);
    
    % project onto positive SI direction
    theta(i) = acos(abs(tmp_principle_axes(SIdirection,:)));
    
end

theta = theta*180/pi;

%% apply smoothing
[vx, vy, vz] = ind2sub(opt.size(1:3), aifROI);
V = [vx, vy, vz];
V = V(:,opt.SIdirection);

for v = unique(V)'
    ind = find(V == v);
    theta(ind) = mean(theta(ind));
end    

end

