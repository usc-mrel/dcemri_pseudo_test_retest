function WC = CELLgetthresh(WC,TH)
%   Extracts wavelet coefficients from a structure and places then in a
%   vector
%   
%   Author: R Marc Lebel
%   Date:   10/2011
%   
%   Usage: wc = CELLgetthresh(SP,SPth,th)
%   
%   Input:
%   WC: (cell) array
%   TH: threshold value cell array
%   
%   Output:
%   WC: column vector

%   Check inputs
if nargin ~= 3
    error('Funtion requires 3 inputs');
end


%   If input is array
if ~iscell(WC) || ~iscell(TH)
    error('Inputs must be cell arrays');
end


%   Determine transform order and extract the number of frames at each one
N = length(WC);

%   Loop through transform order and the various coefficients
for i = 1:N
    WCt = WC{i};
    THt = TH{i};
    fnames = fieldnames(WCt);
    nf = size(fnames);
    for j = 1:nf
        fname = char(fnames(j));
        if ~strcmp(fname,'odd_dims') && ~strcmp(fname,'size') && ~strcmp(fname,'nd')
            coefs = WCt.(fname);
            thrsh = THt.(fname);
            if ~isempty(coefs) && ~isempty(thrsh)
                
                %   Get indices to zero and direction
                ind0 = abs(coefs) < thrsh;
                phase = coefs./(abs(coefs)+eps);
                
                %   Shrink
                coefs = coefs - thrsh.*phase;
                
                %   Zero out
                coefs(ind0) = complex(eps,eps);
                WCt.(fname) = coefs;
            end
        end
    end
    WC{i} = WCt;
end
