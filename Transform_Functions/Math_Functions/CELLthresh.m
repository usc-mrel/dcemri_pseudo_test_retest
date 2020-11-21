function WC = CELLthresh(WC,th,p)
%   Extracts wavelet coefficients from a structure and places then in a
%   vector
%   
%   Author: R Marc Lebel
%   Date:   10/2011
%   
%   Usage: wc = CELLthresh(WC,th)
%   
%   Input:
%   WC: (cell) array
%   th: threshold value (percentile or absolue)
%   p: use percentile
%   
%   Output:
%   WC: column vector

%   Check inputs
if nargin < 1
    error('Funtion requires 1 or more inputs');
end
if nargin < 2
    th = 0.1;
end
if nargin < 3
    p = 1;
end

%   If input is array
if ~iscell(WC)
    error('Input must be a cell array');
end

if p == 1
    %   Obtain column of coefficients
    wc = CELLextract(WC);
    wc = abs(wc);
    
    %   Determine threshold
    th = quantile(wc,th);
    clear wc
end

%   Determine transform order and extract the number of frames at each one
[N1,N2] = size(WC);

%   Loop through transform order and the various coefficients
for k = 1:N1
for i = 1:N2
    WCt = WC{k,i};
    fnames = fieldnames(WCt);
    nf = size(fnames);
    for j = 1:nf
        fname = char(fnames(j));
        if ~strcmp(fname,'odd_dims') && ~strcmp(fname,'size') && ~strcmp(fname,'nd')
            coefs = WCt.(fname);
            if ~isempty(coefs)
                coefs(abs(coefs) < th) = complex(eps,eps);
                WCt.(fname) = coefs;
            end
        end
    end
    WC{k,i} = WCt;
end
end
