function img = iSS(WCD,opt)
%   Applies invere support transform
%   
%   Author: RML
%   Date: 04/2014
%   
%   Usage: img = iSS(WCD,opt)
%   
%   Input:
%   WCD: image of size RO x PE x NS x NT (x NR)
%   opt: Specifies transform options for recon. Use SPSENSE_optset.m
%   	Must include fields:
%           opt.IREF: reference image of size:
%               RO x PE x NS x 1 (x NR)
%           opt.class: data type
%   
%   Output:
%   WC: sparse wavelet coefficient cell

%   Create output
img = zeros(size(WCD),opt.class);

for i = 1:opt.size(4)
    img(:,:,:,i) = WCD(:,:,:,i) .* opt.SPSUP;
end

