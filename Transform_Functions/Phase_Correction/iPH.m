function img = iPH(img,opt)
%   Inverse phase correction
%   
%   Author:
%   R. Marc Lebel
%   05/2011
%   
%   Usage:
%   img = iPH(img,opt)
%   
%	Input:
%   img: complex valued image
%   opt: options structure with field opt.Ph and opt.phase_corr
%   
%	Output:
%   img: complex valued image

%   Apply phase to the image, if desired
if any(opt.phase_corr) && isfield(opt,'Ph') && ~isempty(opt.Ph)
%     PH = conj(opt.Ph);
%     img = img .* PH;
    img = img ./ opt.Ph;
    
    %   Make real (must retain complex nature for calculations)
%     img = complex(real(img),eps*ones(size(img)));
    img = complex(img) + complex(zeros(size(img)), eps*ones(size(img)));
end
