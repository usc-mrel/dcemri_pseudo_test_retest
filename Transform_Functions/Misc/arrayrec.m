function img = arrayrec(img,fc)
%   Optimal image combination for array coil imaging
%   
%   Usage: img2 = arrayrec(img,f)
%   Author: R Marc Lebel
%       Reference: MRM 47:539-548 (2002)
%   Date: 04/2007
%   
%   Inputs:
%   img: individual array coil complex images np x nv (x ns) x nc.
%   fc:  filter cutoff frequency (def. 0.1)
%   
%   Output:
%   img2: Reconstructed images np x nv (x ns)

%   Check inputs
if nargin < 1
    error('arrayrec: function requires at least one input');
end
if nargin < 2 || isempty(fc)
    fc = 0.1;
end

%   Check image size
% np = size(img);
% if length(np) < 3 || length(np) > 4
%     error('arrayrec: image input must be of size np x nv x (x ns) x nc');
% end
[np,nv,ns,nc] = size(img);
if nc == 1
    img = reshape(img,[np nv 1 ns]);
    [np,nv,ns,nc] = size(img);
end

%   Check image data type
if isa(img,'double')
    dt = 0;
elseif isa(img,'single')
    dt = 1;
else
    error('arrayrec: image must be of type double or single');
end

%   Create low resolution (filtered) image
%   Convert to k-space
imgl = fftshift(fft2(single(img)));
if exist('GERecon','file') == 3
    H = GERecon('Windows.KaiserBessel',np,nv,1,0.3,2);
else
    H = butterworth(ones(np,nv),fc);
end
imgl = imgl.*repmat(H,[1 1 ns nc]);
clear H
imgl = ifft2(ifftshift(imgl));

%   Create reference and sensitivity images
%   Could improve somewhat using polynomial fitting for b
alpha = sqrt(sum(abs(imgl).^2,4));
b = imgl./repmat(alpha,[1 1 1 nc]);
clear alpha imgl

%   Compute optimal reconstruction
img = 1./sqrt(sum(abs(b).^2,4)) .* sum(img.*conj(b),4);
clear b

%   Retain only real component
%   could set negative noise to zero, but sometimes this is needed (like subsequent averaging)
%img = real(img);
% img(img < 0) = 0;

end

function Mh = butterworth(M,kc)

%   Check base input argument
if nargin < 1
    error('butter2 requires at least one input')
end

%   Check k space size
N = size(M);
if length(N) < 2 || length(N) > 3
    error('butter2: wrong data size');
end

%   Assign default cutoff frequency
if nargin < 2
    kc = 0.7;
end

%   Create 1D filter
if N(1) == 1 || N(2) == 1
    k1 = -1:2/(length(M)-1):1;
    B = (1./(1+((k1./kc).^2).^4));
elseif ndims(M) == 2
    %   Create 2D filter
    [k1,k2] = meshgrid(-1:2/(N(2)-1):1,-1:2/(N(1)-1):1);
    L = k1.^2 + k2.^2;
    B = 1./(1+(L/kc.^2).^4);
    clear k1 k2 L
elseif ndims(M) == 3
    [k1,k2,k3] = ndgrid(-1:2/(N(1)-1):1,-1:2/(N(2)-1):1,-1:2/(N(3)-1):1);
    L = k1.^2 + k2.^2 + k3.^2;
    B = 1./(1+(L/kc.^2).^4);
end

%   Apply filter
Mh = M.*B;
end