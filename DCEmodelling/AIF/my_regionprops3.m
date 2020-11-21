function  param = my_regionprops3(img, voxelsize, varargin)
% function  param = my_regionprops3(img, varargin)
%   this is my poor man's implementation of Matlab's regionsprop3 to make
%   the STARDCE package work on Matlab < R2017b
%
%
%  Yannick 2020

if isempty(voxelsize)
    voxelsize = [1 1 1];
end

% VoxelIdxList
param = regionprops(img, 'PixelIdxList');
[param.VoxelIdxList] = param.PixelIdxList;
param = rmfield(param,'PixelIdxList');

% Volume
tmp_volume = num2cell(cellfun(@(x) length(x)*prod(voxelsize), {param.VoxelIdxList}));
[param.Volume] = deal(tmp_volume{:});

[Nx, Ny, Nz] = size(img);


if ismember('Zextend', varargin)
    for nvessel = 1:size(param,1)
        [~, ~, vz] = ind2sub([Nx Ny Nz], param(nvessel).VoxelIdxList);
        param(nvessel).Zextend = max(vz) - min(vz);
    end
    
end

if ismember('EigenValues', varargin) || ismember('EigenVectors', varargin)
    for nvessel = 1:size(param,1)
        [vx, vy, vz] = ind2sub([Nx Ny Nz], param(nvessel).VoxelIdxList);
        V = [vx*voxelsize(1), vy*voxelsize(2), vz*voxelsize(3)];
        V = V - repmat(mean(V, 1), [size(V,1) 1]);
        S = (V'*V) / param(nvessel).Volume;
        [V, E] = eig(S);

        [param(nvessel).EigenValues, tmp_I] = sort(diag(E), 'descend');
        param(nvessel).EigenVectors = V(:,tmp_I);
    end
    
end

if ismember('BoundingBox', varargin)
    for nvessel = 1:size(param,1)
        [vx, vy, vz] = ind2sub([Nx Ny Nz], param(nvessel).VoxelIdxList);
        param(nvessel).BoundingBox = [min(vx), max(vx); min(vy), max(vy); min(vz), max(vz)];
    end
    
end

if length(param) == 1
    param = struct2table(param, 'AsArray', true);
else
    param = struct2table(param);
end

end

