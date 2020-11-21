function img = iREF(imgD,opt)
%   Applies the inverse reference image subtraction
%   
%   Author: RML
%   Date: 05/2012
%   
%   Usage: img = iREF(imgD,opt)
%   
%   Input:
%   imgD: difference image of size RO x PE x NS (x NT x NR)
%   opt: Specifies transform options for recon. Use SPSENSE_optset.m
%   	Must include fields:
%           opt.IREF: reference image of size:
%               RO x PE x NS x 1 (x NR)
%           opt.class: data type
%   
%   Output:
%   imgD: difference image (same size)

% %   Create output
% img = zeros(size(imgD),opt.class);
% 
% for i = 1:opt.size(4)
%     img(:,:,:,i,:) = imgD(:,:,:,i,:) + opt.IREF(:,:,:,1,:);
% end

img = imgD;