function FD = fFD(img,FDorder)
%   Forward finite differeneces
%   
%   Author:
%   R Marc Lebel
%   06/2010
%   
%   Usage:
%   FD = fFD(img,opt)
%   
%	Input:
%   img: input image
%   FDorder: cell arry specifying the transform direction
%       ie: {1,2,3}
%   
%	Output:
%   FD: finite difference structure

%   Check inputs
if nargin < 2
    error('fFD requires two inputs');
end

%   Get image size
% [np nv ns nt nr] = size(img);
FD = cell([1,length(FDorder)]);

%   Loop through transform directions
for j = 1:length(FDorder)
    
    if iscell(FDorder)
        drctn = FDorder{j};
    else
        drctn = FDorder(j);
    end
    
    if ~any(drctn == [1 2 3 4])
        error('Unsupported dimension');
    end
    
    %   Perform computation
    FD{j}.fd = fFDMEX(img,drctn);
    
%     FD{j}.X = zeros([np nv ns nt nr]);
%     FD{j}.X(1,:,:,:,:) = img(1,:,:,:,:);
%     for i = 2:np
%         FD{j}.X(i,:,:,:,:) = img(i,:,:,:,:) - img(i-1,:,:,:,:);
%     end
    
end

end
