function [ Ct ] = model_standard(N, Nt, vp, Ktrans, Kep, AIF, t, method)
%function [ Ct ] = model_standard(N, Nt, vp, Ktrans, Kep, AIF, t, method)
%MODEL_STANDARD
%   implements the standard model used for DCE-MRI.
%
% Input:
%   (Nt: number of time points, N: number of voxels)
%   vp      volume of plasma compartment, dimension (N,1)
%   Ktrans  volume transfer constant in [sec^(-1)], dimension (N,1)
%   Kep     = Ktrans/ve, dimension (N,1)
%   AIF     vector of function values (Nt,1) or function handle to AIF
%   t       time vector in sec, dimension (Nt,1)
%   method  integration method
%               'sum': 
%                   piece-wise constant integration / summation of values
%               'trapz':
%                   trapedezoidal integration
%               'adaptive':
%                   adaptive integration (AIF has to be a function handle)
%               'fft':
%                   convolution is performed by applying FFT
%           
% Output:
%   Ct  tracer concentration in tissue, dimension (Nt,N)
%
% References:
%   [1] S.P. Sourbron, D.L. Buckley, "Classic models for dynamic contrast-enhanced MRI",
%       NMR in Biomedicine (2013).
%
%   Yannick Bliesener 2016
%

% Ct = vp*AIF + Ktrans * int_0^t AIF(tau)*exp(-Kep*(t - tau)) dtau

if isa(AIF, 'function_handle')
    AIFval = AIF( t );
else
    if length(t) ~= length(AIF)
        error('model_etofts:Dimension mismatch for t and AIF!');
    end

    AIFval = AIF;
end

vp = reshape(vp, N, 1);
Ktrans = reshape(Ktrans, N, 1);
Kep = reshape(Kep, N, 1);

AIFval = reshape(AIFval, 1, Nt);

t = reshape(t, [1 Nt]);

% check if time bins are uniform, the last one may contain rounding errors
uniform_time = (length(uniquetol(t(2:end-1) - t(1:end-2),0.001)) == 1);
deltat = max(t(2:end-1) - t(1:end-2));

% make sure its causal
t = t - t(1);

Ct = vp*AIFval;

if N == 0
   return 
end

switch method
    case 'sum'
        if uniform_time
            F = exp( repmat(-Kep(:), [1 Nt]) .* repmat(t, [N 1]));
            Ct2 = conv2(1, AIFval, F);
            Ct2 = repmat( (deltat.*Ktrans(:)), [1 Nt]) .* Ct2(:,1:Nt);
            Ct = Ct + Ct2;
        else
            % This is super strange but at least it works
            t = reshape(double(t), [1 Nt]);
            AIFval = reshape(double(AIFval), [1, Nt]);
            vp = reshape(double(vp), 1, N);
            Ktrans = reshape(double(Ktrans), 1, N);
            Kep = reshape(double(Kep), 1, N);
            
            Ct = model_standard_mex(t, AIFval, vp, Ktrans, Kep);
        end
        
    case 'trapz'
        error('MODEL_STANDARD: trapz not maintained')
        
        deltat = t(end) - t(end-1);
        
        parfor k=1:N
            for tt=2:Nt
                % tt = 1 the integral is zeros anyways
                
                F = exp( -Kep(k).*(t(tt)-t(1:tt)));
                Ct(k,tt) =  Ct(k,tt) + (deltat.*Ktrans(k)) * trapz( AIFval(1:tt).*F );
            end
        end
    case 'adaptive'
        error('MODEL_STANDARD: adaptive not maintained')
        
        if ~isa(AIF, 'function_handle')
            error('MODEL_STANDARD: AIF must be a function handle to use adaptive integration!')
        end
        
        parfor k=1:N
            conv_function = @(u,t) Ktrans(k)*exp(-Kep(k)*(t-u)) .* AIF(u);
            integral_function = @(t) integral(@(u) conv_function(u,t), 0, t);
            Ct(:,k) = Ct(:,k) + arrayfun(integral_function,t);
        end
    case 'fft'
        error('MODEL_STANDARD: fft not maintained')
        
        T2 = t( t > 0 );
        Ntt = length(T2);
        
        deltat = t(end) - t(end-1);
        
        [KEP,T] = meshgrid(Kep,T2);
        [KT,~] = meshgrid(Ktrans,T2);
        IRFval = KT.*exp(-T.*KEP);
        clear T KEP KT
        
        AIF_fft = fft( [AIFval; zeros(Ntt-1,1)] );
        IRF_fft = fft( [IRFval; zeros(Nt-1,N)], [], 1 );
        Ct_ifft = ifft( repmat(AIF_fft, [1 N]) .* IRF_fft, [], 1 );
        Ct = Ct + deltat*Ct_ifft(1:Nt,:);
    
    otherwise
        error('MODEL_STANDARD: Unkown method!');
end

Ct = cast(Ct, 'like', AIF);

end

