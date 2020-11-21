function WCD = fSS(img,opt)
%   Applies forward support transform
%   
%   Author: RML
%   Date: 04/2014
%   
%   Usage: imgD = fSS(img,opt)
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
%   WC: sparse wavelet coefficient cell

%   Create output
WCD = zeros(size(img),opt.class);

for i = 1:opt.size(4)
    WCD(:,:,:,i) = img(:,:,:,i) .* opt.SPSUP;
end
