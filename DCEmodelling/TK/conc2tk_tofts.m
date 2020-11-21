function [ Kt, Kep, flag, res, iter ] = conc2tk_tofts(N, Nt, conc, modelParam, opt)
%function [ Kt, Kep, flag, res, iter ] = conc2tk_tofts(N, Nt, conc, modelParam, opt)
%   fits the PK parameters of the Tofts model 
%
%
%   References:
%   1. Murase K. Efficient method for calculating kinetic parameters usingT1-weighted dynamic contrast-enhanced magnetic resonance imaging. Magn. Reson. Med. 2004;51:858?862 doi: 10.1002/mrm.20022.
%
%
%
%   Yannick Bliesener 2017

conc = reshape(conc, [N Nt]);

%% auxiliary variables
init = zeros(2*N,1);
    % initial guess
    
flag_verbose = false;
stepsize     = 1e-3;
maxIter      = 200;

solver       = 'LLSQ';
% solver = 'ncg';

if nargin > 4
    if isfield(opt, 'maxit')
        maxIter = opt.maxit;
    end
    if isfield(opt, 'verbose')
        flag_verbose = opt.verbose;
    end
    if isfield(opt, 'stepsize')
        stepsize = opt.stepsize;
    end
    if isfield(opt, 'init')
        init = opt.init;
    end
    if isfield(opt, 'solver')
        solver = opt.solver;
    end
end


%% Do the deed

switch solver

    case 'ncg'
        
        % Setup inner pcg
        options.method   = 'ncg';
        options.maxIter  = maxIter;
        options.armijoIter = 50;
        options.stepsize = stepsize;
        options.verbose  = false;
        
        [x, flag, output] = fmin(@cost_tofts, @projConstraints, init, options, N, Nt, conc, modelParam, false, [1 1]); 
        
        if flag_verbose
            fprintf('Flag: %i\n', flag)
            fprintf('Residual: %e\n', output.trace(end))
        end
        
        iter    = output.iterations;
        res     = output.trace(end);
%         trace   = output.trace;
        
        Kt  = x(1:N);
        Kep = x(N+1:end);
        
    case 'lm' % levenberg-marquardt
        
        options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt','SpecifyObjectiveGradient',true);

        if flag_verbose
            options.Display   = 'iter';
        else
            options.Display   = 'off';
        end
        
        options.MaxIterations   = maxIter;
        
        [x,~,res,flag,output] = lsqnonlin(@(x) cost_tofts(x,N, Nt, conc, modelParam, true), init, [], [], options);
        
        iter = output.iterations;
        
        if flag_verbose
            fprintf('Report Nelder-Mead\n');
            fprintf('------------------------\n');
            fprintf('Termination flag: %i\n', flag);
            fprintf('Iterations: %i\n', iter);
            fprintf('Function value: %.3e\n', res);
            fprintf('Message: %s\n', output.message);
            fprintf('\n');
        end
        
        Kt  = x(1:N);
        Kep = x(N+1:end);
        
    case 'simplex'
        
        options.MaxIter   = maxIter;
        if flag_verbose
            options.Display   = 'iter';
        else
            options.Display   = 'notify';
        end
        
        [x, f, flag, output] = fminsearch(@(x) cost_tofts(x,N, Nt, conc, modelParam, false), init, options);
        
        iter = output.iterations;
        res = f;
        
        if flag_verbose
            fprintf('Report Nelder-Mead\n');
            fprintf('------------------------\n');
            fprintf('Termination flag: %i\n', flag);
            fprintf('Iterations: %i\n', iter);
            fprintf('Function value: %.3e\n', res);
            fprintf('Message: %s\n', output.message);
            fprintf('\n');
        end
        
        Kt  = x(1:N);
        Kep = x(N+1:end);
        
    case 'LLSQ'
        % Method explained in Ref 1.
        
        Kt  = init(1:N);
        Kep = init(N+1:end);
        
        conc = reshape(conc, N, Nt);
        AIF = reshape(modelParam.AIF, 1, Nt);
        time = reshape(modelParam.time, 1, Nt);
        integration = modelParam.Integration;
        
        intAIF  = integrateAIF(time, AIF, integration);
        intConc = integrateAIF(time, conc, integration);
        
        a = zeros([N Nt 2]);
        a(:, : ,1) = repmat(intAIF, [N 1]);
        a(:, :, 2) = - intConc;
        
        parfor n=1:N
            
            C = reshape(conc(n,:), Nt, 1);
            A = reshape(a(n, :, :), Nt, 2);
            
            if rank(A) < 2
                % essentially Patlak fit w/ vp =0
                A(:,2) = [];
                B = A \ C;
                Kep(n) = 0;
                Kt(n)  = B(1);
            else
                B = A \ C;
                Kep(n) = B(2);
                Kt(n)  = B(1);
            end
        end
        
        flag = 0;
        res = -1;
        iter = -1;
        
    otherwise
        error('Unknown regression method')
end

Kt  = Kt(:);
Kep = Kep(:);

end

function [ cost, grad, H] = cost_tofts(x, N, Nt, conc, modelParam, transposeGrad, gradWeights)

    Kt = x(1:N);
    Kep = x(N+1:end);
    
    Kep( Kep < 0) = 0;

    Ct = model_standard(N, Nt, zeros(size(Kt)), Kt, Kep, modelParam.AIF(:), modelParam.time(:), modelParam.Integration);
    
    deltaC = Ct - conc;
    
    cost = 0.5*norm( deltaC(:), 2)^2;
    
    if nargout > 1
        if nargout > 2
            [ ~, dKtrans, dkep, HI ] = concDeriv_standard(N, Nt, modelParam.time(:), modelParam.AIF(:), Kt, Kep, modelParam.Integration);
        else 
            [ ~, dKtrans, dkep ] = concDeriv_standard(N, Nt, modelParam.time(:), modelParam.AIF(:), Kt, Kep, modelParam.Integration);
        end
        
        grad = zeros(2*N,1);
        grad(1:N) = sum(dKtrans.*deltaC,2);
        grad(N+1:end) = sum(dkep.*deltaC,2);
        
        if nargin > 5
            grad(1:N) = grad(1:N) * gradWeights(1);
            grad(N+1:2*N) = grad(N+1:2*N) * gradWeights(2);
        end
        
        grad = real(grad);
        
        if transposeGrad
            grad = grad';
        end
        
        
        if nargout > 2
           H = zeros(2,2,N);
           
           % dKt^2
           Z = sum( HI.int0.^2 , 2);
           H(1,1,:) = Z;
           % dKep dKt
           Y = repmat(-Kt, [1 Nt]) .* HI.int1;
           Z = sum( Y .* HI.int0, 2) + sum( -HI.int1 .* deltaC, 2);
           H(1,2,:) = Z;
           H(2,1,:) = Z;
           % dKep^2
           Z = sum( Y.^2 , 2 ) +  sum( repmat(Kt, [1 Nt]) .* HI.int2 .* deltaC, 2);
           H(2,2,:) = Z;
        end
        
        
    end
end

function y = projConstraints(x)

    y = x;
    y( x < 0 ) = 0;
%     y( x > 10 ) = 10;

end

