function imgD = fREF(img,opt)
%   Applies the forward reference image subtraction
%   
%   Author: RML
%   Date: 05/2012
%   
%   Usage: imgD = fREF(img,opt)
%   
%   Input:
%   img: image of size RO x PE x NS x NT (x NR)
%   opt: Specifies transform options for recon. Use SPSENSE_optset.m
%   	Must include fields:
%           opt.IREF: reference image of size:
%               RO x PE x NS x 1 (x NR)
%           opt.class: data type
%   
%   Output:
%   imgD: difference image (same size)

%   Create output
imgD = zeros(size(img),opt.class);

if size(opt.IREF,4) == 1
    for i = 1:opt.size(4)
        imgD(:,:,:,i,:) = img(:,:,:,i,:) - opt.IREF(:,:,:,1,:);
    end
else
    imgD = img - opt.IREF;
end

%   Ensure complex
if isreal(imgD)
    imgD(1) = imgD(1) + complex(0,eps);
end
