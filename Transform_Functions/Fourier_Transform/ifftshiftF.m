function f = ifftshiftF(f,dim)
%   Performs an ifftshift in the reciprocal domain
%   
%   Author: Yannick/RML
%   Date: 11/2011
%   
%   Usage: f = ifftshiftF(f,dim)
%   
%   Input:
%   f: input signal
%   dim: specifies the direction to "shift" along
%   
%   Output:
%   F: ifftshifted signal transformed signal

%   Check inputs
if nargin < 1
    error('Function requires at least one input');
end
if nargin < 2
    dim = 1;
end

%   Loop through dimensions
for i = 1:length(dim)
    n = size(f,dim(i));
    ph = exp((-sqrt(-1)*pi*(1:n)));
    rvec = ones(ndims(f),1);rvec(dim(i)) = n;
    ph = reshape(ph,rvec');
    rep = size(f);rep(dim(i)) = 1;
    if isa(f,'single')
        ph = single(ph);
    end
    ph = repmat(ph,rep);
    f = f.*ph;
end
