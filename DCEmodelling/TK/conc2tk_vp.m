function [ Vp ] = conc2tk_vp(N, Nt, conc, AIF, Hct )
% function [ Vp ] = conc2Pk_vp( conc, AIF, Hct )
%CONC2PK_VP 
%   retrieves Vp parameter from concentration and AIF based on vp model
% Input:
%   conc        concentration over time
%                   conc(N, Nt)
%   AIF         arterial input function as vector of values
%   Hct         value of haematocrite
% Output:
%   Vp          fraction plasma volume
%                   Vp
%
%
% Author: Yannick Bliesener 2017

conc = (1 - Hct) * conc;
conc = reshape(conc, [N Nt]);

% discard values of AIF close to zero as these don't allow for PK parameter
% estimation
tmask = ( AIF > eps );

AIF = AIF( tmask );
conc = conc(:, tmask );

Nt = sum( tmask );

if Nt < 2
   error('Insufficient data for Pk parameter estimation') 
end

A = AIF(:);

conc = permute(conc, [2 1]);

x = A \ conc;

Vp = reshape(x(:), [N 1] );

end

