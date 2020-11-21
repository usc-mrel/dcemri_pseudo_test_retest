function [T1,Mo,R2, failureMask] = iDESPOT1(img,opt)
%   Converts images into T1 and Mo maps
%   
%   opt: options structure with fields
%       FA: flip angle (degrees)
%       tr: repetition time (s)
%       B1: transit field scale (fraction)
%       class: specifies 'single' or 'double'

%   Get image size
[np, nv, ns, nFA] = size(img);
fa = pi/180*opt.FA(:);

%   B1
if ~isfield(opt,'B1') || isempty(opt.B1)
    B1 = ones([np, nv, ns],opt.class);
else
    B1 = opt.B1;
end

% if isa(img, 'gpuArray') && existsOnGPU(img)
%     
%     
%     
%     B1 = gpuArray(B1);
%     
%     img = reshape(img, np*nv*ns, nFA);
%     img = transpose(img);
%     
%     b1fa = reshape(fa, nFA, 1) * reshape(B1, 1, np*nv*ns);
%     b = img ./ sin(b1fa);
%     b = reshape(b, [nFA 1 np*nv*ns]);
%     A = img  ./ tan(b1fa);
%     A = cat(2, reshape(A, [nFA 1 np*nv*ns]), ones([nFA 1 np*nv*ns], 'gpuArray'));
%     AH = pagefun(@transpose, A);
%     
%     RHS = pagefun(@times, AH, b);
%     LHS = pagefun(@times, AH, A);
%     AH = repmat(1e-6*eye(2,2, 'gpuArray'), [1 1 np*nv*ns]);
%     LHS = LHS + AH; % just Thikonov regularize, because pinv or their like is not available on GPU....
%     x = pagefun(@inv, LHS, RHS);
%     
%     
%     
%     % very slow
% %     nums = (1:size(A, 3))';
% %     x = arrayfun(@(x)mldivide(A(:,:,x), b(:,:,x)), nums, 'UniformOutput', false);
% %     x = cat(3, x{:});
%     
%     T1 = reshape(x(1,:), [np, nv, ns]);
%     Mo = reshape(x(2,:), [np, nv, ns]);
%     
% else

    %   Precompute some variables
    zna = ones([nFA,1],opt.class);

    %   Initialize output variables
    T1 = zeros([np nv ns],opt.class);
    Mo = zeros([np nv ns],opt.class);
    R2 = zeros([np nv ns],opt.class);

    %   Loop through voxels
    parfor i_np = 1:np
        for i_nv = 1:nv
        for i_ns = 1:ns

            %   Extact signal
            S = img(i_np,i_nv,i_ns,:);
            S = S(:);

            %   Form as matrix problem Ax=b
            b1fa = B1(i_np,i_nv,i_ns).*fa;
            b = S./sin(b1fa);
            A = [S./tan(b1fa) zna];
            
            %   Solve for x
            x = pinv(A) * b;

            %   Store slope and intercept, will convert to T1, Mo after
            T1(i_np,i_nv,i_ns) = x(1);
            Mo(i_np,i_nv,i_ns) = x(2);

        end
        end
    end
% end

%   Regularize solution to prevent divide by zeros
%   Choose the maximum T1 (larger gives less regularization, smaller is more)
T1max = 1000;
mu = 1-exp(-opt.tr/T1max);
T1s = sign(T1-1) .* sqrt((T1-1).^2 + mu^2) + 1;   %   A fix to help avoid divide by zero
ind = T1s<0;

%   Convert to T1
Mo = Mo./(1-T1s);
Mo(ind) = 0;
T1s(ind) = 0;
T1 = -opt.tr./log(T1s);

% failureMask = T1<0 & Mo>0;

%   Fix odd values
ind = ind | isnan(Mo) | isnan(T1) | isinf(Mo) | isinf(T1) | T1<0 | Mo<0;
Mo(ind) = 0;
T1(ind) = 0;

failureMask = ind;

if ~isreal(Mo) || ~isreal(T1)
    warning('Mo or T1 is not real');
    Mo = real(Mo);
    T1 = real(T1);
end

%   Compute R2
if nargout == 3
    imgC = fDESPOT1(T1,Mo,opt);
    SSres = sum((imgC - img).^2,4);
    SStot = sum(img.^2,4);
    R2 = 1 - SSres./SStot;
    R2(R2<0) = 0;
end
