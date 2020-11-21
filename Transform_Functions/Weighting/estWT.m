function WC = estWT(WC,opt)
%   Estimate weighting factors
%   
%   Author:
%   R Marc Lebel
%   11/2011
%   
%   Usage:
%   WT = estWT(WC)
%   
%	Input:
%   WC: coefficients (cell or normal array)
%   opt: options structure with element:
%       opt.l1weight: exponent to the weighting
%       opt.l1reg: regularizer to prevent divide by zeros (1/opt.l1reg is
%                  the maximum possible weight)
%   
%	Output:
%   WT: weights

%   Check inputs
if nargin < 1
    error('Funtion requires one or two inputs');
end

if nargin == 1 || opt.l1weight == 0
    WC = [];
    return
end

%   If input is numeric (non-cell)
if ~iscell(WC)
    WC = 1./(abs(WC).^opt.l1weight + opt.l1reg);
    
elseif iscell(WC)
    
    %   Determine transform order and extract the number of frames at each one
    N = length(WC);
    
    %   Loop through transform order and the various coefficients
    for i = 1:N
        fnames = fieldnames(WC{i});
        nf = size(fnames);
        for j = 1:nf
            fname = char(fnames(j));
            if (~strcmp(fname,'odd_dims') && ~strcmp(fname,'size') && ...
                    ~strcmp(fname,'nd') && ~isempty(WC{i}.(fname)))
                WC{i}.(fname) = 1./(abs(WC{i}.(fname)).^opt.l1weight + opt.l1reg);
            end
        end
    end
    
else
    error('Inconsistent inputs');
end