function WCD = fCD(img,opt)
%   Applies forward baseline subtraction then wavelet
%   
%   Author: RML
%   Date: 05/2012
%   
%   Usage: imgD = fCD(img,opt)
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

if size(opt.IREF,4) == 1
    for i = 1:opt.size(4)
        WCD(:,:,:,i,:) = img(:,:,:,i,:) - opt.IREF(:,:,:,1,:);
    end
else
    WCD = img - opt.IREF;
end

% WCD = fwtN(WCD,opt.wname,opt.worder);
WCD = fVS(WCD,opt);
