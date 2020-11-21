function [ x, exitflag, outputInfo ] = fmin( funObj, funProj, x0, opt, varargin )
%function [ x, exitflag, outputInfo ] = fmin( funObj, projObj, x0, opt, varargin )
%   fmin implements the following algorithms to minimize the function given
%   by funObj:
%       - (projected) gradient descent
%       - (projected) Polyak's Heavyball method
%       - (projected) Nesterov accelerated momentum method
%       - (projected) Adaptive Moment Estimation (Adam)
%       - (projected) Non-linear CG
%
%   Input:
%       funObj     function handle to objective function,
%                   first argument has to be x, others can be passed on
%                   using varargin
%       funProj     function handle to projector onto set of constraints,
%                   argument has to be x
%                   if set to [], this option is disabled
%       x0          initial guess
%       opt         options:
%                       .method     solver:
%                                   - 'graddescent'
%                                   - 'polyak'
%                                   - 'nesterov'
%                       .maxIter    maximum number of iterations
%                       .tol        error tolerance
%                       .errMetric  error metric to measure error tolerance
%                       .armijo     if set to true, Armijo backtracking line search is used to learn stepsize
%                       .stepsize   step size
%                       .verbose    if true, extra information is printed to
%                                   console
%       varargin    arguments passed to funObj
%
%   Output:
%       x           solution vector
%       exitflag    termination flag
%                      -1:  something went wrong
%                       0:  tolerance reached
%                       1:  maximum number of iterations reached
%                       2:  maximum number of line search iteations reached
%                           => could not find descent step
%       outputInfo  Information container
%                       .iterations     number of iterations used
%                       .trace          trace of cost function
%
%   References:
%       [1] Mahdi's lecture (This source code was a homework submission)
%       [2] Kingma, Ba, "ADAM: A METHOD FOR STOCHASTIC OPTIMIZATION", ICLR, 2015.
%       [3] Ruder, "An overview of gradient descent optimization
%       algorithms", arXiv, 2017.
%
% Yannick Bliesener 12/20/16
%
[N,~] = size(x0);

maxIter     = -1;
tol         = -1;
errMetric   = -1;
exitflag    = -1;
epochLength = -1;
    % can be used to learn the stepsize from scratch after 'epochLength'
    % iterations. If epochLength = -1, this option is disabled.
    
LSparam.iter    = 100;
    % Iteration for line search
LSparam.beta    = 0.8;
    % Line search: shrinkage of step size
LSparam.gamma   = 0.5;
    % Line search: desired progress
LSparam.exp     = 1;
    % exponent for inner product of descent direction and gradient

flag_stopCriterion = false;
flag_armijo        = false;
flag_method        = 1;
flag_storeInfo     = false;
flag_projDescent   = false;
flag_verbose       = false;

% initialize a bunch of constants

% poliak & nesterov
eta     = 0;
% Adam
beta    = [0.9 0.999];

% general stuff
stop    = false;
iter    = 0;
x       = x0;

% parse input parameters
if ~isempty(funProj) 
    flag_projDescent = true;
end

if isfield(opt,'maxIter')
    maxIter = opt.maxIter;
    flag_stopCriterion = true;
end

if isfield(opt,'tol') && isfield(opt,'errMetric')
   tol = opt.tol;
   errMetric = opt.errMetric;
   flag_stopCriterion = true;
end

if isfield(opt,'armijo')
    flag_armijo = opt.armijo;
end

if isfield(opt,'armijoIter')
    LSparam.iter = opt.armijoIter;
end

if isfield(opt,'stepsize')
    mu = opt.stepsize;
else
    % no stepsize specified => learn the stepsize
    mu = 1;
    flag_armijo = true;
end

if isfield(opt,'method')
    switch opt.method
        case 'graddescent'
            flag_method = 1;
        case 'polyak'
            flag_method = 2;
            eta = opt.eta;
            xprev = x;
            
%             LSparam.gamma = 0.1;
            
        case 'nesterov'
            flag_method = 3;
            eta = opt.eta;
            xprev = x;
        case 'adam'
            flag_method = 4;
            m = zeros(N,1);
            v = zeros(N,1);
            
            if flag_armijo
                error('Adam optimizer currently does not work with armijo rule...go implement it! ;)')
            end
        case 'ncg'
            flag_method = 5;
            LSparam.gamma = 1e-2;
            LSparam.exp   = 2;
                % I have no idea why but setting this to 2 is way better
                % it's taken from Mikki Lustig's nonlinear cg
                % implementation
            
            flag_armijo = true;
                % ncg always with backtracking line search
            
           	beta = 0;
            v = zeros(N,1);
        otherwise
            error('Unkown optimization algorithm');
    end
end

if isfield(opt,'verbose')
    flag_verbose = opt.verbose;
end

if ~flag_stopCriterion
   error('grad_descent: Please specify stopping criterion!');
end

if nargout > 2
    flag_storeInfo = true;
    trace = [];
end

armijoIter = -1;
LSparam.verbose = flag_verbose;

while ~stop
    
    iter = iter + 1;
   
    switch flag_method
        case 1
            % gradient descent
            [cost, grad] = funObj(x, varargin{:});
            v = -grad;
            
            if flag_armijo
                if mod(iter,epochLength) == 1
                    [mu, armijoIter] = armijo(funObj, funProj, cost, grad, v, x, LSparam, opt.stepsize, varargin{:});
                else   
                    [mu, armijoIter] = armijo(funObj, funProj, cost, grad, v, x, LSparam, mu, varargin{:});
                end
            end
            
            if flag_verbose
                test = zeros(N,1);
            end
            
            % do the deed
            x = x + mu*v;
        case 2
            % polyak
            [cost, grad] = funObj(x, varargin{:});
            v = -grad;
            
            if flag_armijo
                if mod(iter,epochLength) == 1
                    [mu, armijoIter] = armijo(funObj, funProj, cost, grad, v, x, LSparam, opt.stepsize, varargin{:});
                else   
                    [mu, armijoIter] = armijo(funObj, funProj, cost, grad, v, x, LSparam, mu, varargin{:});
                end
                
            end
            
            y = mu*v+eta*(x - xprev);
            
            if flag_verbose
                test = (x - xprev);
            end
            
            % do the deed
            xprev = x;
            x = x + y;
            
        case 3
            % nesterov
            z = x + eta*(x - xprev);
            [cost, grad] = funObj(z, varargin{:});
            v = -grad;
            
            if flag_armijo
                if mod(iter,epochLength) == 1
                    [mu, armijoIter] = armijo(funObj, funProj, cost, grad, v, x, LSparam, opt.stepsize, varargin{:});
                else   
                    [mu, armijoIter] = armijo(funObj, funProj, cost, grad, v, x, LSparam, mu, varargin{:});
                end
            end
            
            if flag_verbose
                test = (x - xprev);
            end
            
            % do the deed
            xprev = x;
            x = z + mu*v;
            
        case 4
            % Adam
            [cost, grad] = funObj(x, varargin{:});
            
            % Mean & Variance
            m = beta(1)*m + (1 - beta(1))*grad;
            v = beta(2)*v + (1 - beta(2))*(grad.^2);
            
            % Bias correction
            mhat = m / (1 - beta(1)^iter);
            vhat = v / (1 - beta(2)^iter);
            
            % do the deed
            xprev = x;
            x = x - mu * mhat ./ (sqrt(vhat) + 1e-8);
        case 5
            % non-linear conjugate gradient
            [cost, grad] = funObj(x, varargin{:});

            if iter > 1
                % FLetcher-Reeves
                beta = (grad'*grad) / (grad_prev'*grad_prev + eps);
            end

            % update search direction
            v = -grad + beta*v;

            % backtracking line search
            if mod(iter,epochLength) == 1
                [mu, armijoIter] = armijo(funObj, funProj, cost, grad, v, x, LSparam, opt.stepsize, varargin{:});
            else
                [mu, armijoIter] = armijo(funObj, funProj, cost, grad, v, x, LSparam, mu, varargin{:});
            end

            % do the deed
            x = x + mu*v;
            grad_prev = grad;

        otherwise
            error('this should not happen!')
    end
    
    if flag_projDescent
       x = funProj(x); 
    end
    
    % check stopping criteria
    if (maxIter >= 0) && (iter >= maxIter)
        stop = true;
        exitflag = 1;
    end
    if (tol > 0) && (errMetric(x) <= tol)
        stop = true;
        exitflag = 0;
    end

    if flag_armijo
        if LSparam.iter == armijoIter
            stop = true;
            exitflag = 2;
        end
    end
    
    if flag_storeInfo
       trace(end+1) = cost; 
    end
    
    if flag_verbose
       fprintf('Iteration: %i\n', iter);
       fprintf('------------------------\n');
       fprintf('Stepsize: %.3e\n', mu);
       if flag_armijo
        fprintf('Armijo-Iterations: %i\n', armijoIter);
       end
       if epochLength > 0
        fprintf('Epoch: %i\n', ceil(iter / epochLength) );
       end
       fprintf('Function value: %.3e\n', cost);
       fprintf('Gradient norm: %.3e\n', norm(v,2));
       
       fprintf('\n');
    end
    
	% adjust step size to make armijo faster
    if flag_armijo
        if armijoIter < 2
            mu = mu / 0.8;
        end
    end
end
    
if flag_storeInfo
    outputInfo.trace = trace;
    outputInfo.iterations = iter;
end

end

function [stepsize, armijoIter] = armijo(funcObj, projObj, curFuncVal, curGrad, curDescentDirection, x, LSparam, startStepsize, varargin)
        
%     normDescentDir = norm( curDescentDirection(:), 2 )^2;
    descentProduct = abs( curGrad(:)'*curDescentDirection(:) )^LSparam.exp;
    stepsize = startStepsize;
    
    if isempty( projObj )
        flag_projDescent = false;
    else
        flag_projDescent = true;
    end
    
    prevPredFuncVal = curFuncVal;
    
    for k=1:LSparam.iter
       
        y = x + stepsize*curDescentDirection;
        if flag_projDescent
            y = projObj(y);
        end
        predFuncVal = funcObj(y, varargin{:});
        
        if LSparam.verbose
            fprintf('Armijo k: %i, stepsize: %e\n',k, stepsize)
            fprintf('Armijo k: %i, current cost: %e\n',k, curFuncVal)
            fprintf('Armijo k: %i, predict cost: %e\n',k, predFuncVal)
        end
        
        if ( predFuncVal < (curFuncVal - LSparam.gamma*stepsize*descentProduct) )
            break;
        elseif (prevPredFuncVal < predFuncVal) && (prevPredFuncVal < curFuncVal)
            stepsize = stepsize / LSparam.beta;
            break;
        end
        
        stepsize = LSparam.beta*stepsize;
        prevPredFuncVal = predFuncVal;
       
    end
    
    armijoIter = k;
end

