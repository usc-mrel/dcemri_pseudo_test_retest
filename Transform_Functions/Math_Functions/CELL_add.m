function WC1 = CELL_add(WC1,WC2)
%   Adds wavelet coefficients
%   
%   Author: R Marc Lebel
%   Date:   01/2011
%   
%   Usage: WC = CELL_add(WC1,WC2)
%   
%   Input:
%   WC1: wavelet coefficient structure
%   WC2: wavelet coefficient structure
%   
%   Output:
%   WC: wavelet coefficient structure, equal to WC1 + WC2

%   Check inputs
if nargin ~= 2
    error('Funtion requires two inputs');
end

%   Determine transform order and extract the number of frames at each one
[N1,N2] = size(WC1);

%   Loop through transform order and the various coefficients
for h = 1:N1
for i = 1:N2
    WCt1 = WC1{h,i};
    WCt2 = WC2{h,i};
    fnames = fieldnames(WCt1);
    nf = size(fnames);
    for j = 1:nf
        fname = char(fnames(j));
        if ~strcmp(fname,'odd_dims') && ~strcmp(fname,'size') && ~strcmp(fname,'nd')
            coefs1 = WCt1.(fname);
            coefs2 = WCt2.(fname);
            if ~isempty(coefs1)
                WCt1.(fname) = coefs1 + coefs2;
            end
        end
    end
    WC1{h,i} = WCt1;
end
end
