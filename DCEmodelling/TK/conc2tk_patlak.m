function [ Kt, Vp ] = conc2tk_patlak(N, Nt, conc, AIF, intAIF, Hct, options )
% function [ Kt, Vp ] = conc2Pk_patlak( conc, time, AIF, intAIF, Hct )
%   retrieves Ktrans and Vp parameter from concentration and AIF/ integral
%   of AIF based on Patlak model
% Input:
%   conc        concentration over time
%                   conc(Nx, Ny, Nz, Nt)
%   AIF         arterial input function as vector of values
%   intAIF      intergal of arterial input function from zero to t
%   Hct         value of haematocrite
% Output:
%   Kt          permeability Ktrans
%                   Kt(Nx, Ny, Nz)
%   Vp          fraction plasma volume
%                   Vp(Nx, Ny, Nz)
%
%
% Author: Yannick Bliesener 07/11/16

flag_enforceConstraints = false;
flag_verbose = false;

if nargin > 6
    if isfield(options,'enforceConstraints')
        flag_enforceConstraints = options.enforceConstraints;
    end
    if isfield(options,'verbose')
        flag_verbose = options.verbose;
    end
end

conc = (1 - Hct) * conc;
conc = reshape(conc, [N Nt]);

AIF = reshape(AIF, [Nt, 1]);
intAIF = reshape(intAIF, [Nt, 1]);

% discard values of AIF close to zero as these don't allow for PK parameter
% estimation
tmask = ( AIF > eps );
Nt = sum( tmask );

AIF = AIF( tmask );
intAIF = intAIF( tmask );

conc = conc(:, tmask );

if Nt < 2
   error('Insufficient data for Pk parameter estimation') 
end

A = [ AIF(:) intAIF(:) ];
conc = permute(conc, [2 1]);

if flag_enforceConstraints
    % solves the fitting as quadratic program to enforce inequality
    % constraints:
    %   Vp >= 0
    %   Kt >= 0
    options =  optimoptions('quadprog','Display','off');
    
    H = A'*A;
    f = -A'*conc;
    
    H = double( H );
    f = double( f );
    
    x = zeros([2 N]);
    for i=1:N
        [x(:,i), fval, exitflag, output] = quadprog(H,f(:,i),[],[],[],[],[0;0],[],[0;0],options);
        if flag_verbose
            fprintf('\tPatlak model solver:\n')
            fprintf('\t\texitflag: %i\n', exitflag)
            fprintf('\t\titerations: %i\n', output.iterations)
        end
    end
else
    % x = pinvA * conc;
    x = A \ conc;
end

Vp = reshape(x(1,:), [], 1);
Kt = reshape(x(2,:), [], 1);

end

