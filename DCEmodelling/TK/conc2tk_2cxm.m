function [ Kt, Vp, Kep, Fp] = conc2tk_2cxm(N, Nt, conc, modelParam, opt)
%function [ Kt, Vp, Kep, Fp] = conc2tk_2cxm(N, Nt, conc, modelParam, opt)
%   fits the PK parameters of the 2-compartment exchange model 
%
%
%   References:
%   1. Flouri D, Lesnic D, Sourbron SP. Fitting the two-compartment model in DCE-MRI by linear inversion. Magn. Reson. Med. 2016;76:998?1006 doi: 10.1002/mrm.25991.
%   2. Sourbron SP, Buckley DL. Classic models for dynamic contrast-enhanced MRI. NMR Biomed. 2013;26:1004?1027 doi: 10.1002/nbm.2940.
%
%   Yannick Bliesener 2019

conc = reshape(conc, [N Nt]);

%% auxiliary variables
init = zeros(4*N,1);
    % initial guess
    
solver       = 'LLSQ';

if nargin > 4
    if isfield(opt, 'solver')
        solver = opt.solver;
    end
end

%% Do the deed

switch solver
        
    case 'LLSQ'
        % Method explained in Ref 1.
        
        Vp  = init(1:N);
        Kt  = init(N+1:2*N);
        Kep = init(2*N+1:3*N);
        Fp  = init(3*N+1:end);
        
        conc = reshape(conc, N, Nt);
        AIF = reshape(modelParam.AIF, 1, Nt);
        time = reshape(modelParam.time, 1, Nt);
        integration = modelParam.Integration;
        
        intAIF  = integrateAIF(time, AIF, integration);
        intConc = integrateAIF(time, conc, integration);
        
        a = zeros([N Nt 4]);
        a(:, : ,1) = - integrateAIF(time, intConc, integration);
        a(:, :, 2) = - intConc;
        a(:, : ,3) = repmat(intAIF, [N 1]);
        a(:, :, 4) = repmat(AIF, [N 1]);
        
        parfor n=1:N
            
            C = reshape(conc(n,:), Nt, 1);
            A = reshape(a(n, :, :), Nt, 4);
            
            if rank(A) == 3
                % essentially Tofts fit, but meaning of variables is different
                % See Ref 2, Table 2
                A(:,[1 2]) = [];
                B = A \ C;
                
                Fp(n)  = Inf;
                Vp(n)  = B(2);
                Kep(n) = 0;
                Kt(n)  = B(1);
            elseif rank(A) == 2
                % essentially vp model fit
                A(:,[1 3]) = [];
                B = A(:,2) \ A(:,1);
                
                Fp(n)  = Inf;
                Vp(n)  = -B;
                Kep(n) = 0;
                Kt(n)  = 0;
            else
                B = A \ C;
                
                T  = B(3) / (B(1)*B(4)); 
                Te = B(2)/B(1) - T;
                Tp = 1 ./ (B(1)*Te);
                
                Fp(n)  = B(4);
                Vp(n)  = B(4)*Tp;
                Ve     = B(4)*(T-Tp);
                PS     = Ve ./ Te;
                Kt(n)  = PS*Fp(n) ./ (PS + Fp(n));
                Kep(n) = Kt(n) ./ Ve;
            end
        end
        
    otherwise
        error('Unknown regression method')
end

Vp  = Vp(:);
Kt  = Kt(:);
Kep = Kep(:);
Fp  = Fp(:);

end