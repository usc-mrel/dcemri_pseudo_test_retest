function wc = CELLextract(WC,skp)
%   Extracts wavelet coefficients from a cell array and places them in a vector
%   
%   Author: R Marc Lebel
%   Date:   10/2011
%   
%   Usage: wc = CELLextract(WC)
%   
%   Input:
%   WC: (cell) array
%   skp: number of entries to skip
%   
%   Output:
%   wc: column vector

%   Check inputs
if nargin < 1
    error('Funtion requires at least 1 input');
end
if nargin < 2
    skp = 1;
end

%   If input is array
if ~iscell(WC)
    wc = WC(1:skp:end);
    
%   If input is cell
else
    
    %   Initialize output
    wc = [];
    
    %   Determine transform order
    [N1,N2] = size(WC);
    
    %   Loop through transform order and the various coefficients
    for h = 1:N1
    for i = 1:N2
        WCt = WC{h,i};
        fnames = fieldnames(WCt);
        nf = size(fnames);
        for j = 1:nf
            fname = char(fnames(j));
            if ~strcmp(fname,'odd_dims') && ~strcmp(fname,'size') && ~strcmp(fname,'nd')
                coefs = WCt.(fname);
                if ~isempty(coefs)
                    coefs = coefs(1:skp:end);
                    wc = [wc;coefs(:)];
                end
            end
        end
    end
    end

end
