function f = iFastFT(f,opt,shift)
%   Computes the [shifted] ifft along the requested dimension(s)
%   
%   Author: RML
%   Date: 09/2011
%   
%   Usage: f = iFastFT(F,opt)
%   
%   Input:
%   F: input signal
%   opt: Specifies transform options. Must contain fields:
%       .dim: specifies transform direction(s)
%       .shift: perform fft shift
%   
%   Output:
%   f: Inverse Fourier transformed signal

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
        f = ifft(ifftshift(f,dim(i)),[],dim(i));
    end

else
    
    %   Loop through dimensions
    for i = 1:length(dim)
        f = ifft(f,[],dim(i));
    end
    
end

% sz = size(f);
% f = f * sqrt(prod(sz(dim)));
