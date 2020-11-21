function [ IC ] = GIC( K, nu, logL )
%function [ IC ] = GIC( K, nu, L )
%   computes the generalized information criterion
%
%   Input:
%       K           number of free parameters
%       logL        maximum value of log-liklihood function

%   Output:
%       IC          guess what that is?!
%
%   References:
%       [1] P. Stoica and Y. Selen, ?Model-order selection,? IEEE Signal Process. Mag., vol. 21, no. 4, pp. 36?47, 2004.
%
% Yannick Bliesener 2017

IC = nu*K - 2*logL;

end
