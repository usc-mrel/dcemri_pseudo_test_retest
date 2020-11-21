function img = iImFilt(imgF,opt)
%   Applies the forward Fourier filter in space
%   
%   Author: RML
%   Date: 10/2014
%   
%   Usage: imgF = iImFilt(img,opt)
%   
%   Input:
%   img: image of size RO x PE x NS x NT (x NR)
%   opt: Specifies transform options for recon. Use SPSENSE_optset.m
%   	 Must include fields:
%        opt.SF_kern: frequency filter of size RO x PE x NS
%   
%   Output:
%   imgF: filtered image

%   Get size
[np, nv, ns, nt, nr] = size(imgF);

%   Take fourier transform of the image
img = fFastFT(imgF,opt.FTdim,1);

%   Multiply with kernel
for t = 1:nt
for r = 1:nr
    img(:,:,:,t,r) = img(:,:,:,t,r) .* opt.SF_kern;
end
end

%   Return to image space
img = iFastFT(img,opt.FTdim,1);

%   Ensure signal remains complex
if ~isreal(imgF) && isreal(img)
    img(1) = complex(img(1),0);
end

