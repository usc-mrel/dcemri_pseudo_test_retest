function Cb = AIF_georgiou(t, param )
% function Cb = AIF_georgiou(t, param )
%AIF_georgiou
% implements a population-based Arterial Input Function
%
% The AIF has been derived in [1]. This is the AIF of a patient with
% Hct=0.35 following a dose = 0.1mmol/kg (patient weight = 72kg), 
% that is closest to the cohort AIF 
%
% Input:
%   t               time vector in sec
%   param           (optional) data container to specify a user defined set of
%                   parameters for the AIF
% Output:
%   Cb              AIF of arterial blood flow (this is blood concentration)
%
% References:
%   1. Georgiou L, Wilson DJ, Sharma N, Perren TJ, Buckley DL. A functional form for a representative individual arterial input function measured from a population using high temporal resolution DCE MRI. Magn Reson Med. 2019;81(3):1955-1963.
%
% Yannick Bliesener 2020
%

% convert to min
t = t / 60;

if ( nargin > 2 ) && ~isempty( param )
    % read user defined input values
    A       = param.A;
    m       = param.m;
    alpha   = param.alpha;
    beta    = param.beta;
    tau     = param.tau;
else
    % default values
    % taken from [1], Table 3
    A     = [0.37, 0.33, 10.06]; % mM
    m     = [0.11, 1.17, 16.02]; % min^{-1}
    alpha = 5.26;
    beta  = 0.032;  % min 
    tau   = 0.129;  % min
end

sizet = size(t);

t = reshape(t, [length(t), 1]);

% Eq. 1 in Ref [1]
Cb = zeros( size(t) );
Cb(t >= 0) = sum( exp(-t(t >= 0) * m) .* A, 2);

n = 0;
while n*tau < t(end)
    indt = (n*tau <= t) & (t < (n+1)*tau);
    gammasum = 0;
    for j = 0:n
        a = (j+1)*alpha+j;
        b = beta;
        x = t(indt) - j*tau;
        gammasum = gammasum + gampdf(x,a+1,b);
    end
    Cb(indt) = Cb(indt) .* gammasum;
    n = n + 1;
end

Cb = reshape(Cb, sizet);

end

