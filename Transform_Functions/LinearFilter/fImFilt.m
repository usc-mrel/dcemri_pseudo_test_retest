function imgF = fImFilt(img,opt)
%   Applies the forward Fourier filter in space
%   
%   Author: RML
%   Date: 10/2014
%   
%   Usage: imgF = fImFilt(img,opt)
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
[np, nv, ns, nt, nr] = size(img);

%   Take fourier transform of the image
imgF = fFastFT(img,opt.FTdim,1);

%   Multiply with kernel
for t = 1:nt
for r = 1:nr
    imgF(:,:,:,t,r) = imgF(:,:,:,t,r) .* opt.SF_kern;
end
end

%   Return to image space
imgF = iFastFT(imgF,opt.FTdim,1);

%   Ensure signal remains complex
if ~isreal(img) && isreal(imgF)
    imgF(1) = complex(imgF(1),0);
end

