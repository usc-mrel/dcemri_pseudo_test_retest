function [ Cb ] = AIF_weinmann( t, groundedDelay, param )
% function [ Cb ] = AIF_weinmann( t, groundedDelay, param )
%AIF_weinmann
% implements a population-based Arterial Input Input function for the generalized pharmacokinetic
% model (e.g. Tofts, Extended-Tofts).
% The AIF has been measured in [1], and fitted in Ref [3].
%
% Input:
%   t               time vector in sec
%   groundDelay     (optional) if set, values at time zero will be grounded, i.e., set
%                   to zero. Default: 0
%   param           (optional) data container to specify a user defined set of
%                   parameters for the AIF
% Output:
%   Cb  AIF of arterial blood flow
%
% References:
%   1. Weinmann HJ, Laniado M, Mützel W. Pharmacokinetics of GdDTPA/dimeglumine after intravenous 
%      injection into healthy volunteers. Physiol Chem Phys Med NMR. 1984;16(2):167-172. http://www.ncbi.nlm.nih.gov/pubmed/6505043.
%   2. Orton MR, Collins DJ, Walker-Samuel S, et al. Bayesian estimation of pharmacokinetic parameters for DCE-MRI with a 
%      robust treatment of enhancement onset time. Phys Med Biol. 2007;52(9):2393-2408.
%   3. Tofts PS, Kermode AG. Measurement of the blood?brain barrier permeability and leakage space using dynamic MR imaging. 1. Fundamental concepts. Magn Reson Med. 1991;17(2):357-367.
%
% Yannick Bliesener 2019
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
    D     = param.D;
    a     = param.a;
    m     = param.m;
else
    % default values
    % taken from [2] which in turn took them from [3]
    D     = 0.1;            % dose 0.1 mmol/Kg
    a     = [3.99 4.78];    % kg/ L
    m     = [0.144 0.0111]; % min^{-1}
end

Cb = zeros( size(t) );

% the following implements Eq. 1.
for i = 1:2
    Cb = Cb + a(i) .* exp(-m(i).*t);
end

Cb = D .* Cb;

if flag_ground
    Cb = Cb .* (t > 0);    
end

end

