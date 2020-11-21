function [ IC ] = BIC( K, N, logL )
%function [ IC ] = BIC( K, N, L )
%   computes the Bayes information criterion
%
%   Input:
%       K           number of free parameters
%       N           number of observations
%       logL        maximum value of log-liklihood function

%   Output:
%       IC          guess what that is?!
%
%   References:
%       [1] K. P. Burnham and D. R. Anderson, ?Multimodel Inference,? Sociol. Methods Res., vol. 33, no. 2, pp. 261?304, 2004.
%       [2] P. Stoica and Y. Selen, ?Model-order selection,? IEEE Signal Process. Mag., vol. 21, no. 4, pp. 36?47, 2004.
%
% Yannick Bliesener 2017

IC = log(N)*K - 2*logL;

end
