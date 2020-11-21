function img = iFD(FD,FDorder)
%   Inverse finite differeneces
%   
%   Author:
%   R Marc Lebel
%   06/2011
%   
%   Usage:
%   img = iFD(FD,opt)
%   
%	Input:
%   FD: gradient cell structure
%   FDorder: cell specifying the transform direction
%       ie: {'1','2','3','4'}

%   Check inputs
if nargin < 2
    error('iFD requires two inputs');
end

%   Initialize output image
img = 0;

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
    img = img + iFDMEX(FD{j}.fd,drctn);
    
%     imgt = zeros(size(FD{j}.X));
%     imgt(1,:,:,:,:) = FD{j}.X(1,:,:,:,:);
%     for i = 2:size(imgt,1)
%         imgt(i,:,:,:,:) = FD{j}.X(i,:,:,:,:) + imgt(i-1,:,:,:,:);
%     end
    
end

end