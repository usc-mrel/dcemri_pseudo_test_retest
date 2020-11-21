function img = iFFTU(k,opt)
%   Computes the [shifted] undersampled ift along the requested dimension(s)
%   
%   Author: RML
%   Date: 09/2011
%   
%   Usage: k = iFFTU(img,opt)
%   
%   Input:
%   k: compressed input signal
%   opt: Specifies transform options. Must contain fields:
%       .dim: specifies transform direction(s)
%       .shift: perform fft shift
%       .size: full size of img
%       .U: sampling indices
%   
%   Output:
%   k: Compacted Fourier transformed signal

%   Check inputs
if nargin < 1
    error('Function requires at two inputs');
end

%   Upsample
img = zeros(opt.size,class(k));
img(opt.U) = k;

%   Perform requested fft
if opt.FTshift
    
    %   Loop through dimensions
    for i = 1:length(opt.FTdim)
        img = ifft(ifftshift(img,opt.FTdim(i)),[],opt.FTdim(i));
    end

else
    
    %   Loop through dimensions
    for i = 1:length(opt.FTdim)
        img = ifft(img,[],opt.FTdim(i));
    end
    
end

% sz = size(img);
% img = img * sqrt(prod(sz(opt.FTdim)));
