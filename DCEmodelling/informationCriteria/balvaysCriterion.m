function [FMI, FRI] = balvaysCriterion(data,estimate)
%function [FMI, FRI] = balvaysCriterion(data,estimate)
%   computes Balvays Criterion for quality of model fit
%
%
%   Output:
%       FMI     fraction of modeled information 
%               (the closer to 1 the better)
%       FRI     fraction of residual information
%               (the closer to zero the better)
%
%   References:
%    [1] D. Balvay et al., "New criteria for assessing fit quality in dynamic 
%        contrast-enhanced T 1-weighted MRI for perfusion and permeability 
%        imaging," Magn. Reson. Med., vol. 54, no. 4, pp. 868?877, 2005.
%
%   Yannick Bliesener 2018

[N, Nt] = size( data );

maxlag = Nt-1;
pOrder = 3;

estimate  = estimate';
data      = data';
residuals = data - estimate;

Pdd0 = Inf(N,1);
Prr0 = Inf(N,1);
for n=1:N
    Rdd = xcorr( data(:,n), maxlag );
    Rdd = Rdd(maxlag+2:end);
    p = polyfit(1:maxlag,Rdd',pOrder);
    Pdd0(n) = real( polyval(p,0) );
    
    Rrr = xcorr( residuals(:,n), maxlag );
	Rrr = Rrr(maxlag+2:end);
    p = polyfit(1:maxlag,Rrr',pOrder);
    Prr0(n) = real( polyval(p,0) );
end

% Eq. 12 in ref 1
errNormSq = Prr0;
d0NormSq  = Pdd0;
resNormSq = sum( abs(residuals).^2, 1 );

% Eq. 13 in ref 1
FMI = 1 - errNormSq(:) ./ d0NormSq(:);
FRI = errNormSq(:) ./ resNormSq(:);

end