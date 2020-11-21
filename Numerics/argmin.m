function [x, f, outputinfo] = argmin(funObj, x0, options, varargin)

if isfield(options, 'MaxIter')
    maxIter = options.MaxIter;
else
    maxIter = 500;
end

if isfield(options, 'PROGTOL')
    tol = options.PROGTOL;
else
    tol = 1e-9;
end

LSparam.iter                        = 100;
LSparam.beta                        = 0.8;
LSparam.gamma                       = 1e-2;
LSparam.exp                         = 2;
LSparam.verbose                     = 0;

x                                   = x0;
t                                   = 1;
f_old                               = 0;
exit                                = false;
opt                                 = varargin{2};
for iter = 1:maxIter
    [f, g] = funObj(x, varargin{:});
    fprintf(['Iteration ' num2str(iter) ', fval: ' num2str(f) '.\n']);
    if iter == 1
        d                           = -g;
        outputinfo.trace            = f;
    else
        gotgo                       = g_old'*g_old;
        beta                        = (g'*g)/(gotgo);
        d                           = -g + beta*d;
        outputinfo.trace(end+1)     = f;
    end
    
    [t, ~]                          = armijo(funObj, f, g, d, x, LSparam, t, varargin{:});
    
    x = x + t*d;
    g_old = g;
    
    if iter == maxIter
        outputinfo.msg              = ['Maximum iteration reached.'];
        exit = true;
    end
    
    if max(abs(t*d)) < tol
        outputinfo.msg              = ['Stepsize below tolerance.'];
        exit = true;
    end
    
    if abs(f_old - f)/f_old < tol && iter > 1
        outputinfo.msg              = ['Function value changing below tolerance.'];
        exit = true;
    end
    
    if exit
        outputinfo.trace(end+1)     = f;
        break;
    end
    
    f_old = f;
end
end

function [stepsize, armijoIter] = armijo(funcObj, curFuncVal, curGrad, curDescentDirection, x, LSparam, startStepsize, varargin)

descentProduct                  = abs( curGrad(:)'*curDescentDirection(:) )^LSparam.exp;
stepsize                        = startStepsize;
prevPredFuncVal                 = curFuncVal;

for k = 1:LSparam.iter
    
    y                           = x + stepsize*curDescentDirection;
    predFuncVal                 = funcObj(y, varargin{:});
    
    if LSparam.verbose
        fprintf('Armijo k: %i, stepsize: %e\n',k, stepsize)
        fprintf('Armijo k: %i, current cost: %e\n',k, curFuncVal)
        fprintf('Armijo k: %i, predict cost: %e\n',k, predFuncVal)
    end
    
    if ( predFuncVal < (curFuncVal - LSparam.gamma*stepsize*descentProduct) )
        break;
    elseif (prevPredFuncVal < predFuncVal) && (prevPredFuncVal < curFuncVal)
        stepsize                = stepsize / LSparam.beta;
        break;
    end
    
    stepsize                    = LSparam.beta*stepsize;
    prevPredFuncVal             = predFuncVal;
end

armijoIter                      = k;
end
