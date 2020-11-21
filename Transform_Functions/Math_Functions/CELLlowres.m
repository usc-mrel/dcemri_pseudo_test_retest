function WC = CELLlowres(WC,n)
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
%   
%   Output:
%   WC: column vector

%   Check inputs
if nargin < 1
    error('Funtion requires 1 or more inputs');
end
if nargin < 2
    n = 1;
end

%   If input is array
if ~iscell(WC)
    error('Input must be a cell array');
end


%   Determine transform order and extract the number of frames at each one
[N1,N2] = size(WC);

%   Loop through transform order and the various coefficients
for k = 1:N1
for i = 1:n
    WCt = WC{k,i};
    fnames = fieldnames(WCt);
    nf = size(fnames);
    for j = 1:nf
        fname = char(fnames(j));
        if ~strcmp(fname,'odd_dims') && ~strcmp(fname,'size') && ~strcmp(fname,'nd')
            coefs = WCt.(fname);
            if ~isempty(coefs)
                coefs = complex(eps,eps)*ones(size(coefs),class(coefs));
                WCt.(fname) = coefs;
            end
        end
    end
    WC{k,i} = WCt;
end
end
