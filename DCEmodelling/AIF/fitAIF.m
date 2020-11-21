function [ param ] = fitAIF( C, t, sigmoidType, initialParam )
%function [ param ] = fitAIF( C, t, sigmoidType, initialParam )
%   fits concentration time curve to AIF function which is either
%   parameterized as:
%       1) a gamma variate plus incomplete gamma sigmoid
%       2) a gamma variate plus sigmoid with ratio of two exponentials
%
%   Input:
%       C               concentration time curve
%       t               time vector
%       sigmoidType     type of sigmoid: 1,2
%       initialParam    container of initial parameters
%                       for incomplete gamma sigmoid (1):
%                           .a      amplitudes [2,1]
%                           .delta 	time delays [2,1]
%                           .alpha 	upslope parameter [2,1]
%                           .tau   	downslope parameter [2,1]
%                           .T      downslope parameter for sigmoid [1,1]
%                       for exponential sigmoid (2):
%                           .a      amplitudes [2,1]
%                           .delta 	time delays [2,1]
%                           .alpha 	upslope parameter [2,1]
%                           .tau   	downslope parameter [2,1]
%                           .beta   sigmoid parameter [1,1]
%   Output:
%       same as initialParam
%
% Yannick 2018
%

verbose = false;

% options = optimoptions('fmincon','Algorithm','trust-region-reflective');
options = optimoptions('fmincon','Algorithm','interior-point');
try
    options.SpecifyObjectiveGradient = true;
catch
    options.GradObj = 'on';
end
options.Display = 'off';

Ng = length(initialParam.a) - 1;
    % number of gamma functions

x0 = struct2x(sigmoidType, initialParam );

% nonnegative constraints
switch sigmoidType
    case 1
        % enforces sigmoid to come after gamma variate 
        A(1,:) = [0 0 1 -1 0 0 0 0 0];
        b(1)   = 0;
        % enforces gamma variate to have higher amplitude than sigmoid
        A(2,:) = [-1 1 0 0 0 0 0 0 0];
        b(2)   = 0;
        % alpha >= tau for gamma variates to give it the characteristic bump
        % shape (see https://en.wikipedia.org/wiki/Gamma_distribution for
        % an illustration of this)
        A(3,:) = [0 0 0 0 -1 0 1 0 0];
        b(3)   = 0;
        % [A delta alpha tau beta]
        lb = [zeros(1,Ng+1) zeros(1,Ng+1) ones(1,Ng) 0 zeros(1,Ng+1) 0];
        ub = [5*ones(1,Ng+1), t(end)*ones(1,Ng+1), Inf(1,Ng+1) Inf(1,Ng+1), Inf];
    case 2
        A = zeros(3*Ng,4*(Ng+1)+1);
        b = zeros(3*Ng,1);
        % enforces sigmoid to come after gamma variate
        for i=1:Ng
            A(i,Ng+1+1)   = 1;
            A(i,Ng+1+1+i) = -1;
        end
        b(1:Ng) = - (t(2)-t(1));
            % make the difference at least one time step to see the
            % difference
        % enforces gamma variate to have higher amplitude than sigmoid
        for i=1:Ng
            A(i+Ng,1)   = -1;
            A(i+Ng,1+i) = 1;
        end
        % alpha >= tau for gamma variates to give it the characteristic bump
        % shape (see https://en.wikipedia.org/wiki/Gamma_distribution for
        % an illustration of this)
%         for i=1:Ng
%             A(i+2*Ng,2*(Ng+1)+i) = -1;
%             A(i+2*Ng,3*(Ng+1)+i) = 1;
%         end
        % beta > alpha for the sigmoid to look like a sigmoid
        A(end+1, :) = [zeros(1,Ng+1) zeros(1,Ng+1) zeros(1,Ng) 1 zeros(1,Ng+1) -1];
        b(end+1) = 0;
        % [A delta alpha tau beta]
        lb = [zeros(1,Ng+1) zeros(1,Ng+1) ones(1,Ng) 0 zeros(1,Ng+1) 0];
        ub = [5*ones(1,Ng+1), t(end)*ones(1,Ng+1), 15*ones(1,Ng+1) Inf(1,Ng+1), 1];
end

inputClass = class(C);
C  = double(C);
t  = double(t);
x0 = double(x0);
A  = double(A);
b  = double(b);
lb = double(lb);
ub = double(ub);

[x,fval,exitflag,output] = fmincon(@(x) AIF( x, C, t, sigmoidType ), x0, A, b, [], [], lb, ub, [], options);

if verbose
    fprintf('Exitflag: %i\n', exitflag)
    fprintf('Iterations: %i\n', output.iterations)
end

switch inputClass
    case 'single'
        x = single(x);
    case 'double'
        x = double(x);
    otherwise
        error('fitAIF: Unkown class')
end

param = x2struct(sigmoidType, x );

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cost, grad] = AIF( x, C, t, sigmoidType )

    param = x2struct( sigmoidType, x);
    
    [ A, deriv ] = AIFgammavariates( t, param );
    
    cost = 0.5*norm( A(:) - C(:) )^2;
    
    if nargout > 1
        
        Ng = (length(x) - 5) / 4;
    
        g = zeros(length(t),4*(Ng+1)+1);
        g(:,1:Ng+1)              = deriv.dA;
        g(:,Ng+2:2*(Ng+1))       = deriv.dDelta;
        g(:,2*(Ng+1)+1:3*(Ng+1)) = deriv.dAlpha;
        g(:,3*(Ng+1)+1:4*(Ng+1)) = deriv.dTau;
        
        switch sigmoidType
            case 1
                g(:,4*(Ng+1)+1) = deriv.dT(:);
            case 2
                g(:,4*(Ng+1)+1) = deriv.dBeta(:);
            otherwise
                error('Unkown sigmoid type')
        end

        grad = g'*(A(:) - C(:));

        
    end
end

function x = struct2x(sigmoidType, param )
    switch sigmoidType
        case 1
            x = [param.a param.delta param.alpha param.tau param.T];
        case 2
            x = [param.a param.delta param.alpha param.tau param.beta];
        otherwise
            error('Unkown sigmoid type')
    end
end

function param = x2struct( sigmoidType, x )

    Ng = (length(x) - 5) / 4;

    switch sigmoidType
        case 1
            param.a     = x(1:Ng+1);
            param.delta = x(Ng+2:2*(Ng+1));
            param.alpha = x(2*(Ng+1)+1:3*(Ng+1));
            param.tau   = x(3*(Ng+1)+1:4*(Ng+1));
            param.T     = x(4*(Ng+1)+1);
        case 2
            param.a     = x(1:Ng+1);
            param.delta = x(Ng+2:2*(Ng+1));
            param.alpha = x(2*(Ng+1)+1:3*(Ng+1));
            param.tau   = x(3*(Ng+1)+1:4*(Ng+1));
            param.beta  = x(4*(Ng+1)+1);
        otherwise
            error('Unkown sigmoid type')
    end

end

