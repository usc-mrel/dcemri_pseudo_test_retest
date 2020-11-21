function imgS = fSENSE(img,opt)
%   Applies the forward SENSE operator to the image
%   
%   Author: RML
%   Date: 11/2011
%   
%   Usage: imgS = fSENSE(img,opt)
%   
%   Input:
%   img: Combined coil image of size RO x PE x NS x NT x 1
%   opt: Specifies transform options for recon. Use SPSENSE_optset.m
%   	Must include fields:
%           opt.S: Sensitivity maps of size:
%               RO x PE x NS x 1 x NR
%   
%   Output:
%   img: image of size RO x PE x NS x NT x NR

% %   Get sizes
% nt = opt.size(4);
% nr = opt.size(5);

%   Make multicoil image
% imgS = repmat(img,[1 1 1 1 nr]).*repmat(opt.S,[1 1 1 nt 1]);

% %   Alternate approach (no repmat)
% imgS = zeros(opt.size,opt.class);
% for i = 1:nt
%     for j = 1:nr
%         imgS(:,:,:,i,j) = img(:,:,:,i).*opt.S(:,:,:,1,j);
%     end
% end

if isa(img, 'gpuArray') && existsOnGPU(img)
    imgS = repmat(img,[1 1 1 1 opt.size(5)]).*repmat(opt.S,[1 1 1 opt.size(4) 1]);
else
    if isfield(opt,'S') && ~isempty(opt.S) && numel(opt.S)~=1
        img = single(img);
        img = complex(img);
        imgS = fSENSE_MEX(img,opt.S,opt.size);
    else
        imgS = img;
    end
end

end
