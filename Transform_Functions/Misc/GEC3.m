function [img,b] = GEC3(img,header,b)
%   GE-like C3 image combination
%   

if nargin < 2
    header.RawHeader.phase_scale = 1;
    header.RawHeader.channel_combine_filter_width = 0.1;
    header.RawHeader.channel_combine_filter_beta = 2;
end

[np,nv,ns,nc] = size(img);
if isfield(header,'RawHeader')
    ps = header.RawHeader.phase_scale;
    fw = header.RawHeader.channel_combine_filter_width;
    fb = header.RawHeader.channel_combine_filter_beta;
else
    ps = header.rdb_hdr_rec.rdb_hdr_phase_scale;
    fw = header.rdb_hdr_rec.rdb_hdr_channel_combine_filter_width;
    fb = header.rdb_hdr_rec.rdb_hdr_channel_combine_filter_beta;
   
end

%   Create low resolution (filtered) image
%   Convert to k-space
imgl = fftshift(fft2(single(img)));
if np > 512 || nv > 512
    H = GERecon('Windows.KaiserBessel',np/2,nv/2,ps,fw,fb);
    H = imresize(H,[np nv]);
else
    H = GERecon('Windows.KaiserBessel',np,nv,ps,fw,fb);
end
%   GE window failure!!
if sum(abs(H(:))) == 0
    H = butterworth(ones(np,nv),fw);
end


imgl = imgl.*repmat(H,[1 1 ns nc]);
clear H
imgl = ifft2(ifftshift(imgl));

%   Create reference and sensitivity images
if nargin < 3 || isempty(b)
    alpha = sqrt(sum(abs(imgl).^2,4))+eps;
    b = imgl./repmat(alpha,[1 1 1 nc]);
end
clear alpha imgl

%   Compute optimal reconstruction
img = sum(img.*conj(b),4);

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
