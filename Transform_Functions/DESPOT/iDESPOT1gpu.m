function [T1,Mo] = iDESPOT1gpu(img,opt)
%   Converts images into T1 and Mo maps
%   
%   opt: options structure with fields
%       FA: flip angle (degrees)
%       tr: repetition time (s)
%       B1: transit field scale (fraction)
%       class: specifies 'single' or 'double'

%   Get image size
img = gpuArray(img);
[np, nv, ns, nFA] = size(img);
fa = gpuArray(pi/180*opt.FA(:));

%   B1
if ~isfield(opt,'B1') || isempty(opt.B1)
    opt.B1 = ones([np, nv, ns],opt.class);
end
B1 = gpuArray(opt.B1);

%   Precompute some variables
zna = ones([nFA,1],opt.class,'gpuArray');

%   Initialize output variables
T1 = zeros([np nv ns],opt.class,'gpuArray');
Mo = zeros([np nv ns],opt.class,'gpuArray');

%   Loop through voxels
for inp = 1:np
for inv = 1:nv
for ins = 1:ns
    
    %   Extact signal
    S = img(inp,inv,ins,:);
    S = S(:);
    
    %   Form as matrix problem Ax=b
    b1fa = B1(inp,inv,ins).*fa;
    b = S./sin(b1fa);
    A = [S./tan(b1fa) zna];
    
    %   Solve for x
    x = pinv(A) * b;
    
    %   Store slope and intercept, will convert to T1, Mo after
    T1(inp,inv,ins) = x(1);
    Mo(inp,inv,ins) = x(2);
    
end
end
end

%   Convert to T1
Mo = gather(Mo./(1-T1));
T1 = gather(-opt.tr./log(T1));

%   Fix abnormal values
ind = isnan(Mo) | isnan(T1) | isinf(Mo) | isinf(T1) | T1<0 | Mo<0;
Mo(ind) = 0;
T1(ind) = 0.001;

