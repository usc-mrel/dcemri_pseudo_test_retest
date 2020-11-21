function [ dvp, dKtrans, dKep, H ] = concDeriv_standard(N, Nt, tModel, AIF, ktMap, kepMap, method)
%function [ dvp, dKtrans, dKep, H ] = concDeriv_standard(N, Nt, tModel, AIF, ktMap, kepMap, method)
%compConcDeriv_standard
% computes the derivative of the concentration w.r.t. vp, Ktrans, and Kep
% for the standard model:
%       C(t) = v_p*AIF(t) + K^{trans}*\int_0^t AIF(\tau)*\exp{-kep*(t-\tau} d\tau
%
% Input:
%   tModel(Nt,1)    time in s (or compatible with fAIF and Ktrans)
%   AIF(Nt,1)       function handle to arterial input function or vector of
%                   function values
%   ktMap(Np)       array of Ktrans values (Np: pixel)
%   kepMap(Np)      array of kep values (Np: pixel)
%
%   method          (optional), default: adaptive
%                   'adaptive':
%                       adaptive integration
%                       AIF needs to be function handle
%                   'sum':
%                       piece-wise constant summation (uniform spacing
%                       only)
%                   'trapz':
%                       trapedezoidal integration (uniform spacing
%                       only)
% Output:
%   dvp(Np,Nt)      derivative of concentration with respect to vp
%                       (Np: pixel, Nt: time)
%   dKtrans(Np,Nt)  derivative of concentration with respect to Ktrans
%                       (Np: pixel, Nt: time)
%   dkep(Np,Nt)      derivative of concentration with respect to kep
%                       (Np: pixel, Nt: time)
%   H               (optional) data container that stores integrals needed
%                   to form Hessians (k=0,1,2):
%                   
%                   H.dims   = [N Nt]
%                   H.int{k} = \int_0^t AIF(\tau) (t-\tau)^k \exp{-(kep)*(t - \tau)} d\tau
%
%
% Yannick Bliesener 05/22/17

%TODO:
% implement fft integration

if nargin > 4
    switch method
        case 'adaptive'
            flag_integrate = 1;
        case 'sum'
            flag_integrate = 2;
        case 'trapz'
            flag_integrate = 3;
        case 'fft'
            flag_integrate = 4;
        otherwise
            error('compSigDiv_standard: Unknown method for integration!')
    end
else
    if isa(AIF, 'function_handle')
        flag_integrate = 1;
    else
        flag_integrate = 3;
    end
end

% auxiliary variables

if isa(AIF, 'function_handle')
    AIFval = AIF( tModel );
else
    if length(tModel) ~= length(AIF)
        error('compSigDiv_standard: Dimension mismatch for t and AIF!');
    end
    if flag_integrate == 1
        error('compSigDiv_standard: AIF needs to be function handle to use adaptive quadrature!');
    end

    AIFval = AIF;
end

switch flag_integrate
    case 1
        deltat = 0; % so the compiler is happy
    case {2,3,4}
        deltat = tModel(end) - tModel(end-1);
end

tModel = reshape(tModel, [1 Nt]);
AIFval = reshape(AIFval, [1 Nt]);

% dvp
dvp = repmat(AIFval, [N 1]);

switch flag_integrate
    
    case 2
        % this is way faster and the integrals do not have to be perfect,
        % it's just for the gradient
        F = exp( repmat(-kepMap(:), [1 Nt]) .* repmat(tModel, [N 1]));
        int0 = deltat * conv2(1, AIFval', F);
        int0 = int0(:,1:Nt);
        int0(:,1) = 0;
        
        F = F .* repmat(tModel, [N 1]);
        delay = sum( tModel == 0 );
        int1 = deltat * conv2(1, AIFval, F);
        int1 = int1(:,delay+1:Nt+delay);
        
%         tic
%         int0 = zeros(N,Nt);
%         int1 = zeros(N,Nt);
% 
%         for tt=2:Nt
%             % tt = 1 the integral is zeros anyways
%             T = repmat(tModel(tt)-tModel(1:tt), [N 1]);
%             
%             F = exp( -repmat(kepMap(:), [1 tt]) .* T );
%             int0(:,tt) = deltat * sum( repmat(AIFval(1:tt), [N 1]) .* F, 2 );
%             int1(:,tt) = deltat * sum( T .* repmat(AIFval(1:tt), [N 1]) .* F, 2 );
%         end
%         toc
        
        if nargout > 3
            
            F = F .* repmat(tModel, [N 1]);
            delay = sum( tModel == 0 );
            int2 = deltat * conv2(1, AIFval, F);
            int2 = int2(:,delay+1:Nt+delay);
            
%             int2 = zeros(N,Nt);
%             for tt=2:Nt
%                 % tt = 1 the integral is zeros anyways
%                 T = repmat(tModel(tt)-tModel(1:tt), [N 1]);
%                 F = exp( -repmat(kepMap(:), [1 tt]) .* T );
%                 int2(:,tt) = deltat * sum( T.^2 .* repmat(AIFval(1:tt), [N 1]) .* F, 2 );
%             end
        end
    case 4

        T2 = tModel( tModel > 0 );
        Ntt = length(T2);
        
        [KEP,T] = meshgrid(kepMap,T2);
        IRFval1 = exp(-T.*KEP);
        IRFval2 = T .* exp(-T.*KEP);
        IRFval3 = T.^2 .* exp(-T.*KEP);
        clear T KEP
        
        AIF_fft   = fft( [AIFval'; zeros(Ntt+1,1)] );
        IRF_fft1  = fft( [IRFval1; zeros(Nt+1,N)], [], 1 );
        IRF_fft2  = fft( [IRFval2; zeros(Nt+1,N)], [], 1 );
        
        int0 = deltat * ifft( repmat(AIF_fft, [1 N]) .* IRF_fft1, [], 1 );
        int1 = deltat * ifft( repmat(AIF_fft, [1 N]) .* IRF_fft2, [], 1 );
        
        int0 = int0(1:Nt,:);
        int0 = int0';
        int1 = int1(1:Nt,:);
        int1 = int1';
        
        int0(:,1) = 0;  % to be the same as the above integrations
        int1(:,1) = 0;
        
        if nargout > 3
           IRF_fft3 = fft( [IRFval3; zeros(Nt+1,N)], [], 1 );
           int2     = deltat * ifft( repmat(AIF_fft, [1 N]) .* IRF_fft3, [], 1 );
           int2     = int2(1:Nt,:);
           int2     = int2';
        end
        
        int2(:,1) = 0; % to be the same as the above integrations
        
    otherwise
        error('compSigDiv_standard: Method for integration not implemented!')
end

dKtrans = int0;
dKep    = -repmat(ktMap(:), [1 Nt]).*int1;

if nargout > 3
   H.dims = [N Nt];
   H.int0 = int0;
   H.int1 = int1;
   H.int2 = int2;
end

if ~( all(isfinite(dvp(:))) && all(isfinite(dKtrans(:))) && all(isfinite(dKep(:))) )
    
    error('compSigDeriv_standard: Computation of derivatives failed! Vp: %i, Kt: %i, Kep: %i', all(isfinite(dvp(:))), all(isfinite(dKtrans(:))), all(isfinite(dKep(:))))
   
end

end
