function [ x, exitflag, outputInfo ]  = fminGaussNewton( funObj, funProj, funHessian, x0, opt, varargin )
%function [ x, exitflag, outputInfo ]  = fminGaussNewton( funObj, funProj, funHessian, x0, opt, varargin )
%   fminGaussNewton implements the following algorithms to minimize the function given
%   by funObj:
%       - (projected) Gauss Newton
%
%   Input:
%       funObj     function handle to objective function,
%                   first argument has to be x, others can be passed on
%                   using varargin
%       funProj     function handle to projector onto set of constraints,
%                   argument has to be x
%                   if set to [], this option is disabled
%       funHessian  function handle to compute new search direction from
%                   Hessian and gradient
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
%       [1] Some calculus class
%
% Yannick 2017
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
flag_storeInfo     = true;
flag_projDescent   = false;
flag_verbose       = false;

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
        case 'gaussnewton'
            flag_method = 1;
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

LSparam.verbose = flag_verbose;

while ~stop
    
    iter = iter + 1;
   
    switch flag_method
        case 1
            % gradient descent
            [cost, grad, H] = funObj(x, varargin{:});
            
            v = -funHessian(H, grad);
            
            if flag_armijo
                if mod(iter,epochLength) == 1
                    [mu, armijoIter] = armijo(funObj, funProj, cost, grad, v, x, LSparam, opt.stepsize, varargin{:});
                else
                    [mu, armijoIter] = armijo(funObj, funProj, cost, grad, v, x, LSparam, mu, varargin{:});
                end
            end
           
            
            % do the deed
            x = x + mu*v;
            
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
    outputInfo.finalStepsize = mu;
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

