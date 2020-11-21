function nrm = l1norm_gpu(WC,WCd,t,sm)
%   Adds two CELL coefficients then takes the l1-norm
%   
%   Author: R Marc Lebel
%   Date:   09/2011
%   
%   Usage: nm = l1norm(WC,WCd,t,sm)
%   
%   Input:
%   WC:  wavelet coefficient structure
%   WCd: wavelet coefficient structure
%   t: scale for WCd
%   sm: smoothing factor for sqrt
%   
%   Output:
%   nrm: norm of WC+t*WCd

%   Check inputs
if nargin < 1
    error('Funtion requires either 1, 3, or 4 inputs');
end
if nargin == 1
    WCd = WC;
    t = 0;
end
if nargin == 2
    error('Function requires either 1, 3, or 4 inputs');
end
if nargin < 4
    sm = 0;
end

%   Make GPU
sm = gpuArray(sm);
t = gpuArray(t);

%   If both inputs are not cells
if ~iscell(WC) && ~iscell(WCd)
%     nrm = l1_normMEX(WC,WCd,t,sm);
    CA_gpu = gpuArray(WC);
    CB_gpu = gpuArray(WCd);
    nrm_gpu = sum(sqrt(abs(CA_gpu(:) + t*CB_gpu(:)).^2 + sm));
    nrm = gather(nrm_gpu);
    
elseif iscell(WC) && iscell(WCd)
    
    %   Determine transform order and extract the number of frames at each one
    [N1,N2] = size(WC);   
    
    %   Loop through transform order and the various coefficients
    nrm = 0;
    for k = 1:N1
    for i = 1:N2
        WCt = WC{k,i};
        WCdt = WCd{k,i};
        fnames = fieldnames(WCt);
        nf = size(fnames);
        for j = 1:nf
            fname = char(fnames(j));
            if (~strcmp(fname,'odd_dims') && ~isempty(WCt.(fname)) && ~isempty(WCdt.(fname)))
%                 nrm = nrm + sum(abs(WCt.(fname)(:) + t*WCdt.(fname)(:)));

                CA_gpu = gpuArray(WCt.(fname));
                CB_gpu = gpuArray(WCdt.(fname));
                nrm_gpu = sum(sqrt(abs(CA_gpu(:) + t*CB_gpu(:)).^2 + sm));
                nrm = nrm + gather(nrm_gpu);
                
            end
        end
    end
    end
    
else
    error('Inconsistent inputs');
end
