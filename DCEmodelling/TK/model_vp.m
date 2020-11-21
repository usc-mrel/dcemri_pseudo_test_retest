function Ct = model_vp(N, Nt, vp, AIF)
%function Ct = model_vp(N, Nt, vp, AIF)
%MODEL_VP
%   implements the vp (only) model used for DCE-MRI.
%
% Input:
%   vp      volume of plasma compartment    
%   AIF     arterial input function
%
% Output:
%   Ct  tracer concentration in tissue
%
% References:
%   [1] S.P. Sourbron, D.L. Buckley, "Classic models for dynamic contrast-enhanced MRI",
%       NMR in Biomedicine (2013).
%
%   Yannick Bliesener 2019

vp = reshape(vp, N, 1);

AIF = reshape(AIF, 1, Nt);

Ct = vp*AIF;

end

