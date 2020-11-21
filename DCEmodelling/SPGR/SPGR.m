function [S, S1, S2, Soff] = SPGR(varargin)
% forward SPGR sequence model
%
%   S    = M0 * \dfrac{sin(alpha) * (1 -  E1)}{1 - E1 * cos(alpha)} * E2 * exp(-i*opt.omega0*deltaB*TE)
%   S1   = \dfrac{sin(alpha) * (1 -  E1)}{1 - E1 * cos(alpha)}
%   S2   = E2
%   Soff = exp(-i*opt.omega0*deltaB*TE)
%
%   examples:
%       1) [S, S1, S2, Soff] = SPGR(M0, R1, opt)
%       2) [S, S1, S2, Soff] = SPGR(M0, R1, R2, opt)
%       3) [S, S1, S2, Soff] = SPGR(M0, R1, R2, deltaB, opt)
%
%   Yannick 2019

M0 = varargin{1};
R1 = varargin{2};
opt = varargin{end};

alpha = opt.FA * pi / 180;

%   B1
if isfield(opt,'B1') && ~isempty(opt.B1)
    alpha = alpha .* opt.B1;
end

S = M0;
E1 = exp(-opt.TR .* R1);

% compatibility with older Matlab versions
rvec  = [ones(1,ndims(E1)-1), size(E1,ndims(E1))];
rvec2 = ones(1,ndims(E1));
sizeE = size(E1);
sizeA = sizeE;
sizeA(1:ndims(alpha)) = size(alpha);
rvec2( sizeA ~= sizeE ) = sizeA(  sizeA ~= sizeE );

if length(alpha) > 1
    S1 = repmat(sin(alpha), rvec) .* repmat((1 - E1), rvec2) ./ (1 - repmat(cos(alpha), rvec) .* repmat(E1, rvec2));
else
    S1 = sin(alpha) .* repmat((1 - E1), rvec2) ./ (1 - cos(alpha) .* repmat(E1, rvec2));
end
if numel(S) == 1
    S = S .* S1;
else
    if isequal(size(S), size(S1))
        S = S .* S1;
    else
        S = repmat(S, rvec) .* S1;
    end
end

S2 = [];
if (nargin > 3)
    R2 = varargin{3};
    if ~isempty(R2)
        S2 = exp(-opt.TE .* R2);
        S = S .* repmat(S2, rvec2);
    end
end

Soff = [];
if (nargin > 4)
    deltaB = varargin{4};
    if ~isempty(deltaB)
        Soff = exp(-sqrt(-1)*opt.omega0*deltaB*opt.TE);
        S = S .* repmat(Soff, rvec2);
    end
end
