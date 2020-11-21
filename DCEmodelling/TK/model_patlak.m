function Ct = model_patlak(N, Nt, vp, kt, AIF, intAIF)
%function Ct = model_patlak(N, Nt, vp, kt, AIF, intAIF)
%MODEL_PATLAK
%   implements the Patlak model used for DCE-MRI.
%
% Input:
%   vp      volume of plasma compartment
%   kt      volume transfer constant in [sec^(-1)]       
%   AIF     arterial input function
%   intAIF  integral of arterial input function
%
% Output:
%   Ct  tracer concentration in tissue
%
% References:
%   [1] S.P. Sourbron, D.L. Buckley, "Classic models for dynamic contrast-enhanced MRI",
%       NMR in Biomedicine (2013).
%
%   Yannick 2019

vp = reshape(vp, N, 1);
kt = reshape(kt, N, 1);

AIF = reshape(AIF, 1, Nt);
intAIF = reshape(intAIF, 1, Nt);

Ct = vp*AIF + kt*intAIF;

end

