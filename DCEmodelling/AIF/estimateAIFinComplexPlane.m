function [AIF, vesselVp] = estimateAIFinComplexPlane(img, opt)
%function AIF = estimateAIFinComplexPlane(img, opt)
%   extracts the AIF from complex data (jointly)
%
%   some adaptation of Ref 1. Although I didn't do the VARPRO stuff....
%
%   see estimateAIF.m
%
%
%   References:
%       1. Klawer EME, van Houdt PJ, Simonis FFJ, et al. Improved repeatability of dynamic contrast-enhanced MRI using the complex MRI signal to derive arterial input functions: a test-retest study in prostate cancer patients. Magn. Reson. Med. 2019:1?12 doi: 10.1002/mrm.27646.
%       2. Simonis FFJ, Sbrizzi A, Beld E, Lagendijk JJW, van den Berg CAT. Improving the arterial input function in dynamic contrast enhanced MRI by fitting the signal in the complex plane. Magn. Reson. Med. 2016;76:1236?1245 doi: 10.1002/mrm.26023.
%       3. Deoni SCL, Peters TM, Rutt BK. Determination of Optimal Angles for Variable Nutation Proton Magnetic Spin-Lattice, T1, and Spin-Spin, T2, Relaxation Times Measurement. Magn. Reson. Med. 2004;51:194?199 doi: 10.1002/mrm.10661.
%
% Yannick Bliesener, bliesene@Usc.edu, 2019

N = prod(opt.size(1:3));
Nt = opt.size(4);

% get initial guess from magnitude data
optI = opt;
% optI.M0 = median(optI.M0);
% optI.R1 = median(optI.R1);
initAIF = iSPGR(abs(img), optI);
indNoError = sum((isnan(initAIF) | isinf(initAIF)), 2) == 0;
initAIF = iFXL(initAIF, optI);
initAIF = initAIF(indNoError(:), :);
if isempty(initAIF)
    AIF = [];
    return
end
initAIF = mean(initAIF, 1);

if sum(indNoError(:)) < 1
    fprintf(opt.fid, 'Could not find valid vessel pixel for complex AIF fit!\n');
    AIF = nan;
    return
end
initAIF = initAIF ./ opt.bloodvp;

% exclude bad voxels
opt.R1 = opt.R1(indNoError(:));
opt.M0 = opt.M0(indNoError(:));
opt.B1 = opt.B1(indNoError(:));
opt.theta = opt.theta(indNoError(:));
img = reshape(img, [N Nt]);
img = img(indNoError(:), :);
opt.size = [size(img,1) 1 1 Nt];

% outlier robust filtering
% opt.M0 = median(opt.M0);
opt.R1 = median(opt.R1);
opt.R2 = median(opt.R2);

% scale
scale = median(opt.M0);
opt.M0 = opt.M0 ./ scale;
img = img ./ scale;

fittingOptions.method   = 'ncg';
fittingOptions.maxIter  = 200;
fittingOptions.stepsize = 1;
fittingOptions.verbose  = false;

if isfield(opt, 'lambda') && ~isempty(opt.lambda) && opt.lambda > 0
    % this approach allows for joint estiamtion of vp and AIF in vessel
    % voxels
    % by anchoring only the mean value of the vp distribution the vessel vp
    % is allowed to follow a Gaussian distribution
    % cost function:
    %   argmin_{vp, AIF} = \|SPGR(vp, AIF) - "img data"\|_2^2 + lambda * \|vp - "nominal blood vp"\|_2^2
    %
   
    weight = ones(size(initAIF));
    
    % set inital values
    vp  = opt.bloodvp*ones(prod(opt.size(1:3)),1);
    AIF = initAIF;
    
    % alternating minimization
    for iter = 1:20
        
        % estimate vessel vp
        [ vp, exitflag, outputInfo ] = fmin( @VPcost, [], vp, fittingOptions, AIF, img, opt);
        
        % update AIF
        opt.vp = vp;
        [ x, exitflag, outputInfo ] = fmin( @AIFcost, [], [1; AIF(:)], fittingOptions, img, opt, weight );
        
        AIF = x(2:Nt+1);
        AIF = reshape(AIF, [1 Nt]);
        AIF = real(AIF);
    end
    
    vesselVp = nan(N,1);
    vesselVp(indNoError(:)) = vp;
else

    % Weighting as introduced in [1,2]. While [2] mentions that the weights of
    % the tail have been reduced until the fit was no longer changing they do
    % not seem to specifiy exactly how this was done
    % this weighting is now introduced as adaptation of [3] (which use a
    % similar idea for relaxometry ...)
    %
    
    weight = max(0.2, initAIF);
    weight = weight ./ norm(weight(:));

    weight = ones(size(weight));
    
    opt.vp = opt.bloodvp;

    [ x, exitflag, outputInfo ] = fmin( @AIFcost, [], [1; initAIF(:)], fittingOptions, img, opt, weight );
    
    AIF = x(2:Nt+1);
    AIF = reshape(AIF, [1 Nt]);
    AIF = real(AIF);
    
    vesselVp = opt.bloodvp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%    Utilities  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [cost, grad] = AIFcost(x, img, opt, weight)
        
        N = prod(opt.size(1:3));
        Nt = opt.size(4);
        
        M0 = x(1);
        conc = x(2:Nt+1);
        conc = reshape(conc, [1 Nt]);
        
        if numel(opt.vp) > 1
            vp = reshape(opt.vp, [N 1]);
            conc = vp * conc;
        else
            conc = opt.vp .* repmat(conc, [N 1]);
        end
        
        [R1, R2] = fFXL(conc, opt);
        if isempty(opt.theta)
            F = 1/3;
        else
            F = reshape(1/3 - 0.5*sin(opt.theta*pi/180).^2, [N 1]);
            F = repmat(F, [1 Nt]); 
        end
        deltaB = opt.TK.chi .* F .* conc;
        
        [S, S1, S2, Soff] = SPGR(opt.M0, R1, R2, deltaB, opt);
        
        % we really need to update the Matlab version on Clipper but then
        % again I also really don't care anymore ...
        if numel(opt.M0) == 1
            rvec = 1;
        else
            rvec = [1 Nt];
        end
        S1 = repmat(opt.M0, rvec) .* S1;
        S  = M0 .* S;
       
        if size(S,1) == 1
            % can happen if no B1 is used
            S = repmat(S, [N 1]);
        end
        
        weight = repmat(weight, [N 1]);

        cost = 0.5*norm(weight(:) .* (S(:) - img(:))).^2;
        
        if nargout > 1
            grad = zeros(size(x));
            
            E1 = exp(-opt.TR.*R1);
            alpha = opt.FA * pi / 180;
            
            %   B1
            if isfield(opt,'B1') && ~isempty(opt.B1)
                alpha = alpha .* opt.B1;
                alpha = reshape(alpha, [N 1]);
            end
            
            if size(S1,1) == 1
                S1 = repmat(S1, [N 1]);
            end

            if size(E1,1) == 1
                E1 = repmat(E1, [N 1]);
            end
%             S2 = repmat(S2, [N 1]);
%             Soff = repmat(Soff, [N 1]);
            
            % grad w.r.t M0
            grad(1) = 0; %sum(sum(conj(S1 .* S2 .* Soff) .* weight.^2 .* (S - img), 2), 1);
            
            % grad w.r.t. C
            grad(2:Nt+1) = ...
                sum(conj(repmat(opt.M0, rvec) .* repmat((sin(alpha)*opt.TR.*( 1 - cos(alpha) ) ), [1 Nt]) .* E1 ./ (1 - E1.*repmat(cos(alpha),[1 Nt])).^2 .* S2 .* Soff .* opt.TK.r1 + ...    % R1 component 
                S1 .* S2 .* Soff .* (-opt.TE*opt.TK.r2) + ...                                                                        % R2 component
                S1 .* S2 .* Soff .* (-sqrt(-1)*opt.omega0*opt.TK.chi*opt.TE*F) ) .* weight.^2 .* (S - img), 1);                                   % off-res component
        
            grad(2:Nt+1) = repmat(conj(M0(:)), [Nt 1]) .* grad(2:Nt+1);
            
            grad(2:Nt+1) = real(grad(2:Nt+1));
            % allow M0 to be complex to account for constant phase
            % offset/shift
        end
    end

    function [cost, grad] = VPcost(vp, AIF, img, opt)
        
        N = prod(opt.size(1:3));
        Nt = opt.size(4);
        
        AIF = reshape(AIF, [1 Nt]);
        vp = reshape(vp, [N 1]);
        
        conc = vp * AIF;
        
        [R1, R2] = fFXL(conc, opt);
        if isempty(opt.theta)
            F = 1/3;
        else
            F = reshape(1/3 - 0.5*sin(opt.theta*pi/180).^2, [N 1]);
            F = repmat(F, [1 Nt]);
        end
        deltaB = opt.TK.chi .* F .* conc;
        
        [S, S1, S2, Soff] = SPGR(opt.M0, R1, R2, deltaB, opt);
        
        % we really need to update the Matlab version on Clipper but then
        % again I also really don't care anymore ...
        if numel(opt.M0) == 1
            rvec = 1;
        else
            rvec = [1 Nt];
        end
        S1 = repmat(opt.M0, rvec) .* S1;
        
        if size(S,1) == 1
            % can happen if no B1 is used
            S = repmat(S, [N 1]);
        end
        
        cost = 0.5*norm((S(:) - img(:)), 2).^2 + 0.5 * opt.lambda * norm(vp(:) - opt.bloodvp, 2).^2;
        
        if nargout > 1            
            E1 = exp(-opt.TR.*R1);
            alpha = opt.FA * pi / 180;
            
            %   B1
            if isfield(opt,'B1') && ~isempty(opt.B1)
                alpha = alpha .* opt.B1;
                alpha = reshape(alpha, [N 1]);
            end
            
            if size(S1,1) == 1
                S1 = repmat(S1, [N 1]);
            end

            if size(E1,1) == 1
                E1 = repmat(E1, [N 1]);
            end
            
            % grad w.r.t. vp
            grad = sum(repmat(AIF, [N, 1]) .* conj(repmat(opt.M0, rvec) .* repmat((sin(alpha)*opt.TR.*( 1 - cos(alpha) ) ), [1 Nt]) .* E1 ./ (1 - E1.*repmat(cos(alpha),[1 Nt])).^2 .* S2 .* Soff .* opt.TK.r1 + ...    % R1 component
                S1 .* S2 .* Soff .* (-opt.TE*opt.TK.r2) + ...                                                                        % R2 component
                S1 .* S2 .* Soff .* (-sqrt(-1)*opt.omega0*opt.TK.chi*opt.TE*F) ) .* (S - img), 2);                                   % off-res component
            
            grad = real(grad);
            
            % penalty term
            grad = grad + opt.lambda * (vp(:) - opt.bloodvp);
        end
    end
end