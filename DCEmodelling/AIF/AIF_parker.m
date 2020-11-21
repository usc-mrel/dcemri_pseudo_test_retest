function [ Cb ] = AIF_parker( t, groundedDelay, param )
% function [ Cb ] = AIF_parker( t )
%AIF_PARKER
% implements a population-based Arterial Input Input function for the generalized pharmacokinetic
% model (e.g. Tofts, Extended-Tofts).
% The AIF has been derived in [1].
%
% Input:
%   t               time vector in sec
%   groundDelay     (optional) if set, values at time zero will be grounded, i.e., set
%                   to zero. Default: 0
%   param           (optional) data container to specify a user defined set of
%                   parameters for the AIF
% Output:
%   Cb  AIF of arterial blood flow (this is blood concentration)
%
% References:
%   [1] G.J.M. Parker, et al., "Experimentally-Derived Functional Form for a Population-Averaged 
%       High-Temporal-Resolution Arterial Input Function for Dynamic
%       Contrast-Enhanced MRI", Magnetic Resonance in Medicine 56:993?1000, 2006.
%
% Yannick Bliesener 2015
%

if nargin < 2
   flag_ground = 1;
else
   flag_ground = groundedDelay;
end

% convert to min
t = t / 60;

if ( nargin > 2 ) && ~isempty( param )
    % read user defined input values
    A     = param.A;
    T     = param.T;
    sigma = param.sigma;
    alpha = param.alpha;
    beta  = param.beta;
    s     = param.s;
    tau   = param.tau;
else
    % default values
    % values are taken from Table 1.
    A     = [0.809 0.330];
    T     = [0.17046 0.365];
    sigma = [0.0563 0.132];
    alpha = 1.050;
    beta  = 0.1685;
    s     = 38.078;
    tau   = 0.483;
end

Cb = zeros( size(t) );

% the following implements Eq. 1.
for i = 1:2
    Cb = Cb + (A(i)/(sigma(i)*sqrt(2*pi))) * exp(-(t-T(i)).^2 / (2*sigma(i)^2) );
end

Cb = Cb + alpha*exp(-beta*t) ./ ( 1 + exp(-s*(t-tau)) );

if flag_ground
    Cb = Cb .* (t > 0);    
end

end

