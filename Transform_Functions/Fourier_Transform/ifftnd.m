function kout = ifftnd(k,fftdim,shiftdim)

% kout = ifftnd(k,fftdim,shiftdim)
%
% modified IFFT version specifically shift designated dimension
% default(if no shiftdim input), shift both domain
% shiftdim=0, perform only fft
% specific shiftdim, shift only afterwards
%
% Yi Guo 04/28/2013
% modified 05/06/2013
% edited by ZYH, 06/05/2014


if ~exist('shiftdim','var')
    for n = 1: length(fftdim)
        k=ifft(k,[],fftdim(n));
    end
end

if exist('shiftdim','var')
    % set shiftdim=0 to perform only fft
    if shiftdim==1
        for n=1:length(fftdim)
        k = ifftshift(ifft(ifftshift(k,fftdim(n)),[],fftdim(n)),fftdim(n));
        end
    else
        for n=1:length(fftdim)
            k = ifft(k,[],fftdim(n));
        end
    end

    % shift specific dim afterwards
    % if shiftdim
    % for n = 1: length(shiftdim)
    %   k = fftshift(ifft(k,[],fftdim(n)),shiftdim(n));
    % end
    % end

end

kout=k;

end
