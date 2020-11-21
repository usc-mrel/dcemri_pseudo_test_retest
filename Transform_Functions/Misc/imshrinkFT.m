function img2 = imshrinkFT(img,sz2,fc)
%   Interpolates images in 3D

%   Check inputs
if nargin < 3
    fc = 0.9;
end
if nargin < 2
    error('Requires at lest two inputs');
end

%   Get image size
[np, nv, ns, nt] = size(img);
sz1 = [np, nv, ns, nt];

%   Check desired interpolation
switch numel(sz2)
    case 1
        sz2 = round(s.*sz1);
        sz2(4) = sz1(4);
    case 3
        sz2(4) = sz1(4);
    case 4
        sz2(4) = sz1(4);
    otherwise
        error('Wrong interpolation size');
end

%   Transform to k-space
img2 = fFastFT(img,[1 2 3],1);

%   Create interpolation points
ind1 = (round(sz1(1)/2) - round(sz2(1)/2) + 1): (round(sz1(1)/2) + floor(sz2(1)/2));
ind2 = (round(sz1(2)/2) - round(sz2(2)/2) + 1): (round(sz1(2)/2) + floor(sz2(2)/2));
ind3 = (round(sz1(3)/2) - round(sz2(3)/2) + 1): (round(sz1(3)/2) + floor(sz2(3)/2));
ind4 = 1:nt;
if ns == 1
    img2 = img2(ind1,ind2,1   ,ind4);
else
    img2 = img2(ind1,ind2,ind3,ind4);
end

%   Create filter
if fc < 1 && fc > 0
    F = fermi(ones(sz2(1:3),class(img)),fc,1/25);
    for i = 1:nt
        img2(:,:,:,i) = img2(:,:,:,i).*F;
    end
end

%   Convert back to image space
img2 = iFastFT(img2,[1 2 3],1);

%   Scale img2 due to FT
img2 = (prod(sz2)/prod(sz1)) * img2;