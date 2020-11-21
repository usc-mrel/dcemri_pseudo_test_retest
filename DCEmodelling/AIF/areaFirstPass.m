function [AUC1, AIFparam] = areaFirstPass(time, AIF, ibolus, method)
% function AUC1 = areaFirstPass(time, AIF, tbolus)
%   compute area under the first pass of the AIF
%
%   Inputs:
%       time    time vector
%       AIF     vascular input function
%       ibolus  time index of bolus arrival
%       method  'numeric' | 'gaussianfit'
%   Outputs:
%       AUC         the number you came for
%       AIFparam    if the estimate of AUC arose due to fitting of a
%                   Gaussian curve to the first pass, then this is the parameter
%                   container
%
%
%
% Yannick Bliesener, bliesene@Usc.edu, 2020

AIFparam = [];

if nargin < 3
    dtime = max(time(2:end-1) - time(1:end-2)); 
    ibolus = estimateBAT([1 1 1 length(AIF)], AIF, dtime, 'LL');
end

if nargin < 4
    method = 'numeric';
end

switch method
    
    case 'numeric'
        [~, Imax] = max(AIF);

        localMin = (AIF(1:end-2) > AIF(2:end-1)) & (AIF(3:end) > AIF(2:end-1));
        localMin = [false localMin false];
        localMin = find( localMin & (1:length(AIF) > Imax) & ((time - time(Imax) <= 15)) ); % find minimum within 10s of the peak

        if isempty(localMin)
            % find any point within 10s second window
            [~, localMin] = min(abs(time(Imax+1:end) - time(Imax) - 15));
            localMin = localMin + Imax;
        end
        localMin = localMin(1);

        integration_interval = ibolus:localMin;
        AUC1 = trapz(time(integration_interval), AIF(integration_interval));
    case 'gaussianfit'
        
        fitOpt = struct();
        fitOpt.time = time - time(1);
        fitOpt.AIF.paramterization = '2Gammas+Sigmoid';
        fitOpt.AIF.initialParameters = [];
        AIFparam = parameterizeAIF( AIF, fitOpt );
        
        AIFparam.a(2:end) = 0;
        
        normFactor = ( exp(1)./(AIFparam.tau(1)*(AIFparam.alpha(1)-1)) ).^(AIFparam.alpha(1)-1);
        gammaFactor = 1 / (gamma(AIFparam.alpha(1))*AIFparam.tau(1)^AIFparam.alpha(1));
        AUC1 = AIFparam.a(1) * normFactor / gammaFactor; % mM * s = mmol s/L

    otherwise
        error('Unkown method!')
end
end

