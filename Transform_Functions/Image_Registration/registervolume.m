function [imgR,tform,Mot] = registervolume(img,img0,FOV,tform,type,mode, init_tform)

if nargin < 7 || isempty(mode)
    init_tform = affine3d();
end
if nargin < 6 || isempty(mode)
    mode = 'multimodal';
end
if nargin < 5 || isempty(type)
    type = 'rigid';
end
if nargin < 4 || isempty(tform)
    tform = [];
end
if nargin < 3 || isempty(FOV)
    FOV = [1 1 1];
end

interpolator = 'nearest';    % 'nearest', 'linear', 'cubic'

%   Matrix size
[np, nv, ns] = size(img);
dx = FOV(1)/np;
dy = FOV(2)/nv;
dz = FOV(3)/ns;
ref = imref3d([np nv ns],dx,dy,dz);


if isempty(tform)
    %   Estimate
    %   Default parameters
    [optimizer, metric] = imregconfig(mode);
    if strcmp(mode,'multimodal')
        metric.UseAllPixels = 0;
        metric.NumberOfHistogramBins = 512;
        %metric.NumberOfSpatialSamples = 300000;
        metric.NumberOfSpatialSamples = ceil(np*nv*ns/15);
        optimizer.InitialRadius = 1e-5;
        optimizer.GrowthFactor = 1.05;
        optimizer.Epsilon = 1e-6;
    end
    optimizer.MaximumIterations = 1024;
    tform = imregtform(img, ref, img0, ref, type, optimizer, metric,...
        'DisplayOptimization',false,...
        'PyramidLevels',3,...
        'InitialTransformation',init_tform);
    Mot = [tform.T(4,1);tform.T(4,2);tform.T(4,3);...
        180/pi*tform.T(2,1);180/pi*tform.T(3,1);180/pi*tform.T(3,2)];
end


%   Apply
if isreal(img)
    imgR = imwarp(img,ref,tform,interpolator,'OutputView',ref);
else
    imgR = imwarp(real(img),ref,tform,interpolator,'OutputView',ref);
    imgI = imwarp(imag(img),ref,tform,interpolator,'OutputView',ref);
    imgR = complex(imgR, imgI);
end

end

