function img = fPH(img,opt)
%   Forward phase correction
%   
%   Author:
%   R Marc Lebel
%   05/2011
%   
%   Usage:
%   img = fPH(img,opt)
%   
%	Input:
%   img: real valued image
%   opt: options structure with field opt.Ph and opt.phase_corr
%   
%	Output:
%   img: complex valued image

%   Apply phase to the image, if desired
if any(opt.phase_corr) && isfield(opt,'Ph') && ~isempty(opt.Ph)
    img = opt.Ph .* img;
end
