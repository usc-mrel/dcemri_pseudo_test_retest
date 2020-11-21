function PH = estPH(img,opt)
%   Estimate phase correction factor
%   
%   Author:
%   R Marc Lebel
%   08/2011
%   
%   Usage:
%   PH = estPH(img,opt)
%   
%	Input:
%   img: complex valued image
%   opt: options structure
%   
%	Output:
%   PH: phase image

%   Get image size
[np, nv, ns, nt, nr] = size(img);

%   Define filter properties
if ~isfield(opt,'Ph_cf')
    opt.Ph_cf = [0.25 0.25 0.25];
end
T = 1/200;

%   Initialize low-pass Fermi filter
F = ones([np,nv,ns,nt,nr],opt.class);

for i = 1:length(opt.phase_corr)

    switch opt.phase_corr(i)
        
        %   First dimension
        case 1
            if np > 1
                x = -1:2/(np-1):1;
            else
                x = 0;
            end
            F1 = 1./(1+exp((abs(x)-opt.Ph_cf(1))/T));
            F = F .* repmat(F1(:),[1 nv ns nt nr]);
        
        %   Second dimension
        case 2
            if nv > 1
                x = -1:2/(nv-1):1;
            else
                x = 0;
            end
            F1 = 1./(1+exp((abs(x)-opt.Ph_cf(2))/T));
            F = F .* repmat(reshape(F1,[1 nv]),[np 1 ns nt nr]);
        
        %   Third dimension
        case 3
            if ns > 1
                x = -1:2/(ns-1):1;
            else
                x = 0;
            end
            F1 = 1./(1+exp((abs(x)-opt.Ph_cf(3))/T));
            F = F .* repmat(reshape(F1,[1 1 ns]),[np nv 1 nt nr]);
            
        otherwise
            error('Unsupported direction');
    end
end

%   FFT shift the filter
for i = 1:length(opt.phase_corr)
    F = fftshift(F,opt.phase_corr(i));
end

%   Covert image to k-space
PH = fFastFT(img,opt.phase_corr,0);
    
%   Apply filter and extract phase, convert to complex exponential
PH = PH .* F; clear F;
PH = iFastFT(PH,opt.phase_corr,0);
PH = exp(sqrt(-1)*angle(PH));
