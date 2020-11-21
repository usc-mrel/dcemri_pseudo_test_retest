function [ IC ] = BIC_gaussian( K, N, RSS )
%function [ IC ] = BIC_gaussian( K, N, RSS )
%   computes the Bayes information criterion for the special case of
%   Gaussian noise and LS fit
%
%   Input:
%       K           number of free parameters
%       N           number of observations
%       RSS         Residual sum of squares (final value of LS)

%   Output:
%       IC          guess what that is?!
%
%   References:
%       [1] K. P. Burnham and D. R. Anderson, ?Multimodel Inference,? Sociol. Methods Res., vol. 33, no. 2, pp. 261?304, 2004.
%       [2] P. Stoica and Y. Selen, ?Model-order selection,? IEEE Signal Process. Mag., vol. 21, no. 4, pp. 36?47, 2004.
%
% Yannick Bliesener 2017

IC = log(N)*K + N*log( RSS / N + eps );

end
