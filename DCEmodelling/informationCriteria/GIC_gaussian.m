function [ IC ] = GIC_gaussian( K, N, nu, RSS )
%function [ IC ] = BIC_gaussian( K, N, RSS )
%   computes the generalized information criterion for the special case of
%   Gaussian noise and LS fit
%
%   Input:
%       K           number of free parameters
%       N           number of observations
%       nu          penalty factor (the higher the less overfittig)
%       RSS         Residual sum of squares (final value of LS)

%   Output:
%       IC          guess what that is?!
%
%   References:
%       [1] P. Stoica and Y. Selen, ?Model-order selection,? IEEE Signal Process. Mag., vol. 21, no. 4, pp. 36?47, 2004.
%
% Yannick Bliesener 2017

IC = nu*K + N*log( RSS / N + eps );

end

