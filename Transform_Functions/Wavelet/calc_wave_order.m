function worder = calc_wave_order(img,minsize)
%   Estimate wavelet transform order
%   
%   Author:
%   R Marc Lebel
%   07/2011
%   
%   Usage:
%   worder = calc_wave_order(img,minsize)
%   
%	Input:
%   img: N-D image
%   minsize (optional): minimum coefficient length (def: 5)
%   
%	Output:
%   worder: wavelet transform dimensions for use with fwtN.m and iwtN.m


%   Check input arguments
if nargin < 2
    minsize = 5;
end

%   Get image size
sz = size(img);

%   Estimate transform size
ts = floor(log2(sz) - log2(minsize));
mt = max(ts);

%   Check if image is too small
if all(ts < 1)
    warning('calc_wave_order:Image_too_small','Image is too small to be transformed');
    worder = [];
    return;
end

%   Loop through transforms
for i = 0:mt-1
    co = mt - i;
    flag = ts >= co;
    dimstr = [];
    for j = 1:length(sz)
        if flag(j)
            dimstr = [dimstr j]; %#ok<*AGROW>
        end
    end
    worder{i+1} = dimstr;
end

%   Flip string
worder = fliplr(worder);
