function [ IC ] = AIC_gaussian( K, N, RSS, corrected )
%function [ IC ] = AIC_gaussian( K, N, RSS, corrected )
%   computes the Akaike information criterion for the special case of
%   Gaussian noise and LS fit
%
%   Input:
%       K           number of free parameters
%       N           number of observations
%       RSS         Residual sum of squares (final value of LS)
%       corrected   (optional)
%                       0: (default) no correction for finite sample size
%                       1: correction for finite sample size
%   Output:
%       IC          guess what that is?!
%
%   References:
%       [1] K. P. Burnham and D. R. Anderson, ?Multimodel Inference,? Sociol. Methods Res., vol. 33, no. 2, pp. 261?304, 2004.
%       [2] H. Akaike, ?A New Look at the Statistical Model Identification,? IEEE Trans. Automat. Contr., vol. 19, no. 6, pp. 716?723, 1974.
%       [3] P. Stoica and Y. Selen, ?Model-order selection,? IEEE Signal Process. Mag., vol. 21, no. 4, pp. 36?47, 2004.
%
% Yannick Bliesener 2017

if nargin == 3
    corrected = 0;
end

IC = 2*K + N*log( RSS / N + eps );

if corrected
    IC = IC + 2*K*(K+1) / (N - K - 1);
end

end

