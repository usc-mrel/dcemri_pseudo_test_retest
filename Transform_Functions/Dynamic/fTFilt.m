function imgF = fTFilt(img,opt)
%   Applies the forward Fourier filter along the time dimension
%   
%   Author: RML
%   Date: 07/2014
%   
%   Usage: imgF = fTFilt(img,opt)
%   
%   Input:
%   img: image of size RO x PE x NS x NT (x NR)
%   opt: Specifies transform options for recon. Use SPSENSE_optset.m
%   	 Must include fields:
%        opt.TF_kern: frequency filter of size 1 x 1 x 1 x nt x 1
%   
%   Output:
%   imgF: filtered image

%   Get size
[np, nv, ns, nt, nr] = size(img);

%   Check filter
if ~isfield(opt,'TF_kern') || numel(opt.TF_kern) ~= nt
    error('Missing or bad filter');
end

%   Reshape kernel, convert class, replicate to full image size
kern = reshape(opt.TF_kern,[1 1 1 nt 1]);
if strcmp(opt.class,'single') && ~isa(kern,'single')
    kern = single(kern);
elseif strcmp(opt.class,'double') && ~isa(kern,'double')
    kern = double(kern);
end
kern = repmat(kern,[np nv ns 1 nr]);

%   Take fourier transform of the image along time
%   Do not worry about shifting - this means that DC will be at the edge
%   and the filter needs to be designed appropriately
imgF = fft(img,[],4);

%   Multiply with kernel
imgF = imgF .* kern;
clear kern

%   Return to image space
imgF = ifft(imgF,[],4);

%   Ensure signal remains complex
if ~isreal(img) && isreal(imgF)
    imgF(1) = complex(imgF(1),0);
end
