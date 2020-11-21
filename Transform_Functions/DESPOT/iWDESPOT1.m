function [T1,Mo,R2, failureMask] = iWDESPOT1(img,opt)
%   Converts images into T1 and Mo maps
%
%   uses the weighted cost function as described in 
%   ?1. Deoni SCL, Peters TM, Rutt BK. Determination of Optimal Angles for Variable Nutation Proton Magnetic Spin-Lattice, T1, and Spin-Spin, T2, Relaxation Times Measurement. Magn. Reson. Med. 2004;51:194?199 doi: 10.1002/mrm.10661.
%   
%   opt: options structure with fields
%       FA: flip angle (degrees)
%       tr: repetition time (s)
%       B1: transit field scale (fraction)
%       class: specifies 'single' or 'double'
%
%   Yannick 2020

%   Get image size
[np, nv, ns, nFA] = size(img);
fa = pi/180*opt.FA(:);

%   B1
if ~isfield(opt,'B1') || isempty(opt.B1)
    B1 = ones([np, nv, ns],opt.class);
else
    B1 = opt.B1;
end

%   Precompute some variables
zna = ones([nFA,1],opt.class);

%   Initialize output variables
T1 = zeros([np nv ns],opt.class);
Mo = zeros([np nv ns],opt.class);
R2 = zeros([np nv ns],opt.class);

%   Loop through voxels
for i_np = 1:np
    for i_nv = 1:nv
    for i_ns = 1:ns

        %   Extact signal
        S = img(i_np,i_nv,i_ns,:);
        S = S(:);

        %   Form as matrix problem Ax=b
        b1fa = B1(i_np,i_nv,i_ns) .* fa;
        
        b = S./sin(b1fa);
        A = [S./tan(b1fa) zna];
        
        %   Solve for x
        x = pinv(A) * b;

        %   weighting
        W = S ./ max(S);

        fittingOptions.method   = 'ncg';
        fittingOptions.maxIter  = 30;
        fittingOptions.stepsize = 1e-2;
        fittingOptions.verbose  = false;

        [ x, exitflag, outputInfo ] = fmin( @WDESPOT1cost, [], [x(2) / (1-x(1)); x(1)], fittingOptions, S, W, b1fa);

        %   Store slope and intercept, will convert to T1, Mo after
        T1(i_np,i_nv,i_ns) = x(2);              % E1
        Mo(i_np,i_nv,i_ns) = x(1) * (1- x(2));  % M0 * (1 - E1)

    end
    end
end

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cost, grad] = WDESPOT1cost(x, img, W, b1fa)
        
        M0 = x(1);
        E1 = x(2);
        
        S = M0 .* (1-E1) .* sin(b1fa) ./ (1 - E1 .* cos(b1fa));
        
        diff = W(:).*(S(:) - img(:));
        cost = 0.5*norm(diff).^2;
        
        if nargout > 1
            grad = zeros(size(x));
            
            % grad w.r.t M0
            grad(1) = sum(W(:) .* ((1-E1) .* sin(b1fa) ./ (1 - E1 .* cos(b1fa))) .* diff);
            
            % grad w.r.t. E1
            grad(2) = sum(W(:) .* (cos(b1fa) - 1) ./ (1 - E1 .* cos(b1fa)).^2 .* diff);
            
            grad(2) = real(grad(2));
            % allow M0 to be complex to account for constant phase
            % offset/shift
        end
    end

end