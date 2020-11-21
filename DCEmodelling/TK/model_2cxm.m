function [ Ct ] = model_2cxm(N, Nt, vp, ve, Ktrans, Fp, AIF, t, method)
%function [ Ct ] = model_2cxm(N, Nt, vp, ve, Ktrans, Fp, AIF, t, method)
%MODEL_2CXM
%   implements the two-compartment exchange model used for DCE-MRI.
%
% Input:
%   (Nt: number of time points, N: number of voxels)
%   vp      volume of plasma compartment, dimension (N,1)
%   ve      volume of EES compartment, dimension (N,1)
%   Ktrans  volume transfer constant in [sec^(-1)], dimension (N,1)
%   Fp    	CBF, dimension (N,1)
%   AIF     vector of function values (Nt,1) or function handle to AIF
%   t       time vector in sec, dimension (Nt,1)
%   method  integration method
%               'sum': 
%                   piece-wise constant integration / summation of values
%               'trapz':
%                   trapedezoidal integration
%           
% Output:
%   Ct  tracer concentration in tissue, dimension (Nt,N)
%
% References:
%   [1] S.P. Sourbron, D.L. Buckley, "Classic models for dynamic contrast-enhanced MRI",
%       NMR in Biomedicine (2013).
%   [2] Barnes SR, Ng TSC, Santa-Maria N, Montagne A, Zlokovic B V, Jacobs RE. ROCKETSHIP: a flexible and modular software tool for the planning, processing and analysis of dynamic MRI studies. BMC Med. Imaging 2015;15:19.
%
%  Yannick Bliesener 2016

if isa(AIF, 'function_handle')
    AIFval = AIF( t );
else
    if length(t) ~= length(AIF)
        error('model_etofts:Dimension mismatch for t and AIF!');
    end

    AIFval = AIF;
end

vp = reshape(vp, N, 1);
ve = reshape(ve, N, 1);
Ktrans = reshape(Ktrans, N, 1);
Fp = reshape(Fp, N, 1);


AIFval = reshape(AIFval, 1, Nt);

t = reshape(t, [1 Nt]);

% check if time bins are uniform, the last one may contain rounding errors
uniform_time = (length(uniquetol(t(2:end-1) - t(1:end-2),0.001)) == 1);
deltat = max(t(2:end-1) - t(1:end-2));

% make sure its causal
t = t - t(1);

if Ktrans>= Fp
    PS = 10^8; %large value, prevents Inf and NaN errors
else
    PS = Ktrans.*Fp./(Fp-Ktrans);
end
E = PS./(PS+Fp);
e = ve./(vp+ve);
tau_plus  = (E-E.*e+e)./(2.*E).*(1+sqrt(1-(4.*E.*e.*(1-E).*(1-e))./(E-E.*e+e).^2));
tau_minus = (E-E.*e+e)./(2*E).*(1-sqrt(1-(4.*E.*e.*(1-E).*(1-e))./(E-E.*e+e).^2));
k_plus  = Fp./((vp+ve).*tau_minus);
k_minus = Fp./((vp+ve).*tau_plus);
F_plus  =  1*Fp.*(tau_plus-1)./(tau_plus-tau_minus);
F_minus = -1*Fp.*(tau_minus-1)./(tau_plus-tau_minus);


switch method
    case 'sum'
        
        if uniform_time
            F = F_plus.*exp( repmat(-k_plus(:), [1 Nt]) .* repmat(t, [N 1])) + F_minus*exp(repmat(-k_minus(:), [1 Nt]) .* repmat(t, [N 1]));
            Ct = conv2(1, AIFval, F);
            Ct = deltat * Ct(:,1:Nt);
        else
            % variable frame rate
            % This is super strange but at least it works
            t = reshape(double(t), [1 Nt]);
            AIFval = reshape(double(AIFval), [1, Nt]);
            vp = reshape(double(vp), 1, N);
            ve = reshape(double(ve), 1, N);
            Ktrans = reshape(double(Ktrans), 1, N);
            Fp = reshape(double(Fp), 1, N);
            
            error('model_2cxm_mex: Needs to be debugged as this will likely crash!')
            Ct = model_2cxm_mex(t, AIFval, vp, ve, Ktrans, Fp);
        end

        
    case 'trapz'
        % Pre-alocate for speed
        Ct = zeros(N,numel(T1));
        
        parfor k=1:N
            for tt=2:Nt
                
                % The time for T
                T = T1(1:tt);
                CP= Cp(1:tt);
                
                F = CP.*(F_plus*exp(-k_plus*(T(end)-T)) + F_minus*exp(-k_minus*(T(end)-T)));
                
                if(numel(T) == 1)
                    %need this as trapz interprets non array as
                    %Y,DIM input instead of X,Y
                    Ct(k,tt) = 0;
                else
                    Ct(k,tt) = trapz(T,F);
                end
                
                if isnan(Ct(k,tt))
                    Ct(k,tt) = 0;
                end
            end
        end
    
    otherwise
        error('MODEL_2CXM: Unkown method!');
end

Ct = cast(Ct, 'like', AIF);


end

