function WC = absgradMEX(WC,smooth)
%   MEX file to compute the approximate gradient of all wavelet coefficients
%   
%   Author: R Marc Lebel
%   Contact: mlebel@gmail.com
%   Date:   11/2010
%   
%   Usage: WC2 = absgradMEX(WC,smooth)
%   
%   Input:
%   WC: numeric array
%   smooth: small smoothing factor to prevent Inf
%   
%   Output:
%   WC2: numeric array