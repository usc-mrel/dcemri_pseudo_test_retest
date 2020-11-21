function f = fFastFT(f,opt,shift)
%   Computes the [shifted] fft along the requested dimension(s)
%   
%   Author: RML
%   Date: 09/2011
%   
%   Usage: F = fFastFT(f,opt)
%   
%   Input:
%   f: input signal
%   opt: Specifies transform options. Must contain fields:
%       .dim: specifies transform direction(s)
%       .shift: perform fft shift
%   
%   Output:
%   F: Fourier transformed signal

%   Check inputs
if nargin < 1
    error('Function requires at least one input');
end
if nargin < 2
    dim = 1;
    shift = 1;
end
if nargin == 2 && isstruct(opt)
    if ~isfield(opt,'FTshift')
        opt.FTshift = 1;
    end
    if ~isfield(opt,'FTdim')
        opt.FTdim = 1;
    end
    shift = opt.FTshift;
    dim = opt.FTdim;
elseif nargin >= 2
    dim = opt;
end


%   Perform requested fft
if shift
    
    %   Loop through dimensions
    for i = 1:length(dim)
        f = fftshift(fft(f,[],dim(i)),dim(i));
    end

else
    
    %   Loop through dimensions
    for i = 1:length(dim)
        f = fft(f,[],dim(i));
    end
    
end

% sz = size(f);
% f = f / sqrt(prod(sz(dim)));
