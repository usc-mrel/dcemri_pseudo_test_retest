function kVS = fVS(k,opt)
%   Applies the forward view sharing kernel to the k-space
%   
%   Author: RML
%   Date: 08/2011
%   
%   Usage: kG = fVS(k,opt)
%   
%   Input:
%   k: k-space of size RO x PE x NS x NT x NR
%   opt: Specifies transform options for recon. Use SPSENSE_optset.m
%   	Must include fields:
%           opt.vs_kern: view_sharing kernel of size:
%               1 x 1 x 1 x nt x 1
%   
%   Output:
%   kVS: view shared k-space

%   Get size
[np nv ns nt nr] = size(k);

%   Subtract input signal by modifying the central kernel pixel
% lk = ceil(numel(opt.vs_kern)/2);
% opt.vs_kern(lk) = opt.vs_kern(lk)-1;

%   Reshape the kernel and convert class
kern = reshape(opt.vs_kern(1:end),[1 1 1 length(opt.vs_kern) 1]);
if strcmp(opt.class,'single') && ~isa(kern,'single')
    kern = single(kern);
elseif strcmp(opt.class,'double') && ~isa(kern,'double')
    kern = double(kern);
end

%   Replicate k on edges
lk = floor(numel(kern)/2);
kVS = cat(4,repmat(k(:,:,:,1,:),[1 1 1 lk 1]),k,repmat(k(:,:,:,nt,:),[1 1 1 lk 1]));

%   Convolve kernel with k-space
kVS = convnc(kVS,kern,'valid');
% kVS = convnc(k,kern,'same');

%   Ensure signal remains complex
if ~isreal(k) && isreal(kVS)
    kVS(1) = complex(kVS(1),eps);
end
