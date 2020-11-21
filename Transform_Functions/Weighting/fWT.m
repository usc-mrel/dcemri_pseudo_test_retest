function WC = fWT(WCwt,WC)
%   Apply weighting factor
%   
%   Author:
%   R Marc Lebel
%   11/2011
%   
%   Usage:
%   WC = fWT(WCwt,WC)
%   
%	Input:
%   WCwt: weights (cell or normal array)
%   WC: coefficients to be weighted
%   
%	Output:
%   WC: weighted coefficients

if isempty(WCwt) || any(size(WCwt) ~= size(WC))
    return;
end

if iscell(WCwt) && iscell(WC)
    WC = CELL_mult(WC,WCwt);
else
    WC = WCwt.*WC;
end
