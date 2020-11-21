function k = fFFTU(img,opt)
%   Computes the [shifted] undersampled fft along the requested dimension(s)
%   
%   Author: RML
%   Date: 09/2011
%   
%   Usage: k = fFFTU(img,opt)
%   
%   Input:
%   img: fully sampled input signal
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
    error('Function requires at least one input');
end
if nargin < 2
    opt.FTdim = 1;
    opt.FTshift = 1;
    opt.U = 1:numel(img);
    opt.size = size(img);
end

%   Perform requested fft
if opt.FTshift
    
    %   Loop through dimensions
    for i = 1:length(opt.FTdim)
        img = fftshift(fft(img,[],opt.FTdim(i)),opt.FTdim(i));
    end

else
    
    %   Loop through dimensions
    for i = 1:length(opt.FTdim)
        img = fft(img,[],opt.FTdim(i));
    end
    
end

% sz = size(img);
% img = img / sqrt(prod(sz(opt.FTdim)));

%   Undersample
k = img(opt.U);
