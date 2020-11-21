function WC = absgrad_gpu(WC,smooth)
%   Computes the approximate gradient of all wavelet coefficients
%   
%   Author: R Marc Lebel
%   Contact: mlebel@gmail.com
%   Date:   11/2010
%   
%   Usage: WC2 = absgrad(WC,smooth)
%   
%   Input:
%   WC: (cell) array
%   smooth: small smoothing factor to prevent Inf
%   
%   Output:
%   WC2: (cell) array

%   Check inputs
if nargin ~= 2
    error('Funtion requires 2 inputs');
end

smooth = gpuArray(smooth);

%   If input is array
if ~iscell(WC)
    WC_gpu = gpuArray(WC);
    WC_gpu = WC_gpu.*(conj(WC_gpu).*WC_gpu + smooth).^(-0.5);
    WC = gather(WC_gpu);
    
%   If input is cell
else
    
    %   Determine transform order
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
                    WC_gpu = gpuArray(coefs);
                    WC_gpu = WC_gpu.*(conj(WC_gpu).*WC_gpu + smooth).^(-0.5);
                    coefs = gather(WC_gpu);
                    WCt.(fname) = coefs;
                end
            end
        end
        WC{k,i} = WCt;
    end
    end

end
