function sig = fSPGR(varargin)
% forward SPGR sequence model
%
%   examples:
%       1) sig = fSPGR(R1, opt)
%       2) sig = fSPGR(R1, R2, opt)
%       3) sig = fSPGR(R1, R2, deltaB, opt)
%
%   Yannick 2019

R1 = varargin{1};
opt = varargin{end};

alpha = opt.FA * pi / 180;

%   B1
if isfield(opt,'B1') && ~isempty(opt.B1)
    alpha = alpha .* opt.B1;
end

rvec = [ones(1,ndims(R1)-1), size(R1,ndims(R1))];
    % compatibility with older Matlab versions

sig = repmat(opt.M0 .* sin(alpha), rvec)  .* (1 - exp(-opt.TR*R1)) ./ (1 - repmat(cos(alpha),rvec).*exp(-opt.TR*R1));

ind = repmat((cos(alpha) == 1), rvec) | (R1 == 0);
sig(ind) = 0;

if (nargin > 2)
    R2 = varargin{2};
    if ~isempty(R2)
        sig = sig .* exp(-opt.TE*R2);
    end
end

if (nargin > 3)
    deltaB = varargin{3};
    if ~isempty(deltaB)
        sig = sig .* exp(-1i*opt.omega0*opt.TE*deltaB);
    end
end

end

