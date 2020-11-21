function [ C, deriv ] = AIFgammavariates( t, param )
%function [ C ] = AIFgammavariates( t, param )
% computes the AIF based on gamma variates plus a sigmoid and its
% derivative with respect to the parameters
%
%   Input
%       t       time vector in seconds
%       param   (optional) parameter container (see below)
%                   if not specified default values are taken from [1]
%
%   Output
%       C       vector with AIF
%       deriv   container with derivatives 
%
%   References
%       [1] J. U. Fluckiger, M. C. Schabel, and E. V. R. DiBella, ?Model-based blind estimation of
%           kinetic parameters in Dynamic Contrast Enhanced (DCE)-MRI,? 
%           Magn. Reson. Med., vol. 62, no. 6, pp. 1477?1486, 2009.
%       [2] G. J. M. Parker et al., ?Experimentally-derived functional form for a population-averaged 
%           high-temporal-resolution arterial input function for dynamic contrast-enhanced MRI,? 
%           Magn. Reson. Med., vol. 56, no. 5, pp. 993?1000, 2006.
%
%   Author
%       Yannick Bliesener
%   Date
%       03/02/2017
%

flag_computeDerivatives = 0;

sigmoid = false;
sigmoidType = 1;

if nargout > 1
    flag_computeDerivatives = 1;
end
 
if (nargin > 1) && ~isempty( param )
   a     = param.a;
   alpha = param.alpha;
   tau   = param.tau;
   delta = param.delta;
   
   if length(a) > 1
       sigmoid = true;
       if isfield(param, 'T')
           sigmoidType = 1;
           T = param.T;
       elseif isfield(param, 'beta')
           sigmoidType = 2;
           beta = param.beta;
       else
           error('Cannot determine sigmoid type!')
       end
   end
else
   sigmoid = true;
   a     = [6 1.1208 0.3024 0.7164];
   alpha = [2.92 2.92 2.92 2.92];
   tau   = [0.0442 0.143 0.143 0.143] * 60;
   delta = [1 1.2227 1.6083 1.6083] * 60 - 60;
   T     = 7.8940 * 60;
end

C = zeros( size(t) );

if sigmoid 
    Ng = length( a ) - 1;
else
    Ng = 1;
end

%% gamma variates

for i=1:Ng
    tt = t - delta(i);
    
    stepFunc = (tt >= 0);
    tt = tt .* stepFunc;
    
    % auxiliary variables
    normFactor = ( exp(1)./(tau(i)*(alpha(i)-1)) ).^(alpha(i)-1);
    linearTerm = tt.^(alpha(i)-1);
    expTerm    = exp(-tt./tau(i));
    
    % function evalutation
    gammaFunc = normFactor .* linearTerm .* expTerm;
    gammaFunc( ~stepFunc ) = 0;
    C = C + a(i) .* gammaFunc;
    
    % derivatives
    if flag_computeDerivatives
        
        % amplitude
        dA = gammaFunc;
       
        % delta
        dDelta = normFactor .* linearTerm .* expTerm;
        dDelta = dDelta .* ( (1-alpha(i))./(tt + eps) + 1 / tau(i) );
        dDelta( ~stepFunc ) = 0;
        
        % alpha
        dAlpha = normFactor .* linearTerm .* expTerm;
        dAlpha = dAlpha .* (log( (tt + eps) .* exp(1) ./ tau(i) ./ (alpha(i)-1) ) - 1);
        dAlpha( ~stepFunc ) = 0;
        
        % tau
        dTau   = normFactor .* linearTerm .* expTerm;
        dTau   = dTau .* ( (1-alpha(i))./tau(i) + tt ./ tau(i)^2 );
        dTau( ~stepFunc ) = 0;
        
        %
        deriv.dA(:,i)     = dA;
        deriv.dDelta(:,i) = a(i) .* dDelta;
        deriv.dAlpha(:,i) = a(i) .* dAlpha;
        deriv.dTau(:,i)   = a(i) .* dTau;
    end
end


%% sigmoid
if sigmoid 

    switch sigmoidType

        case 1
            % Sigmoid taken from ref [1]

            % auxiliary variables
            tt         = t - delta(Ng+1);
            stepFunc   = (tt >= 0);
            tt         = tt .* stepFunc;

            scaleTime  = (1/tau(Ng+1) - 1/T);
            scaleTime  = max(scaleTime,0);
            normFactor = T / (T - tau(Ng+1)) * (tau(Ng+1)/T)^(-tau(Ng+1)/(T - tau(Ng+1)) );
            expTerm    = exp(-tt./T);
            gammaTerm  = gammainc( scaleTime * tt, alpha(Ng+1), 'lower');

            % function evalutation
            sigmoid  = normFactor .* expTerm .* gammaTerm;
            sigmoid( ~stepFunc ) = 0;
            C = C + a(Ng+1) .* sigmoid;

            % derivatives
            if flag_computeDerivatives

                % amplitude
                dA = sigmoid;

                % delta
                dGamma_dDelta = -scaleTime .* ( scaleTime .* tt ).^(alpha(Ng+1) - 1) .* exp( -scaleTime .* tt );
                dDelta = normFactor  .* expTerm .* ( gammaTerm ./ T + dGamma_dDelta / gamma( alpha(Ng+1) ) );
                dDelta( ~stepFunc ) = 0;

                % alpha
                dGamma_dAlpha = derivGammaInc( alpha(Ng+1), scaleTime.*tt );
                dAlpha = normFactor .* expTerm .* ( -gammaTerm .* psi( alpha(Ng+1) ) + dGamma_dAlpha / gamma(alpha(Ng+1)) );
                dAlpha( ~stepFunc ) = 0;

                % tau
                dGamma_dTau = ( scaleTime .* tt ).^(alpha(Ng+1) - 1) .* exp( -scaleTime .* tt ) .* (-tt./tau(Ng+1).^2);
                dTau = normFactor .* expTerm .* (gammaTerm ./ (T - tau(Ng+1)) + ...
                    gammaTerm .* (-T./(T - tau(Ng+1)).^2 .* log(tau(Ng+1)/T) - 1./(T - tau(Ng+1)) ) + ...
                    dGamma_dTau / gamma(alpha(Ng+1)) );
                dTau( ~stepFunc ) = 0;

                % T
                dGamma_dTau = ( scaleTime .* tt ).^(alpha(Ng+1) - 1) .* exp( -scaleTime .* tt ) .* (tt./T.^2);
                dT = normFactor .* expTerm .* (gammaTerm .* -tau(Ng+1) ./ (T - tau(Ng+1)) ./ T + ...
                    gammaTerm .* ( tau(Ng+1)/(T-tau(Ng+1)).^2 .* log(tau(Ng+1)/T) + tau(Ng+1)/(T-tau(Ng+1))/T) + ...
                    gammaTerm .* tt ./ T^2 + ...
                    dGamma_dTau / gamma(alpha(Ng+1)) );
                dT( ~stepFunc ) = 0;

                %
                deriv.dA(:,Ng+1)     = dA;
                deriv.dDelta(:,Ng+1) = a(Ng+1) .* dDelta;
                deriv.dAlpha(:,Ng+1) = a(Ng+1) .* dAlpha;
                deriv.dTau(:,Ng+1)   = a(Ng+1) .* dTau;
                deriv.dT             = a(Ng+1) .* dT;

            end

        case 2

            % sigmoid taken from ref [2]

            %% Sigmoid
            % auxiliary variables
            tt         = t - delta(Ng+1);
            stepFunc   = (tt >= 0);
            tt         = tt .* stepFunc;

            expTermAlpha = exp(-alpha(Ng+1) .* tt);
            expTermBeta  = exp(-beta .* (tt - tau(Ng+1)) );
            sigmoid  = expTermAlpha ./ ( 1 + expTermBeta );
            sigmoid( ~stepFunc ) = 0;
            C = C + a(Ng+1) .* sigmoid;

            % derivatives
            if flag_computeDerivatives

                % amplitude
                dA = sigmoid;

                % delta
                dDelta = alpha(Ng+1) .* expTermAlpha .* ( 1 + expTermBeta ) - beta .* expTermAlpha .* expTermBeta;
                dDelta = dDelta ./ ( 1 + expTermBeta ).^2;
                dDelta( ~stepFunc ) = 0;

                % alpha
                dAlpha = -tt .* sigmoid;
                dAlpha( ~stepFunc ) = 0;

                % beta
                dBeta = expTermAlpha .* expTermBeta .* (tt - tau(Ng+1)) ./ ( 1 + expTermBeta ).^2;
                dBeta( ~stepFunc ) = 0;

                %tau
                dTau = expTermAlpha .* expTermBeta .* -beta ./ ( 1 + expTermBeta ).^2;
                dTau( ~stepFunc ) = 0;

                %
                deriv.dA(:,Ng+1)     = dA;
                deriv.dDelta(:,Ng+1) = a(Ng+1) .* dDelta;
                deriv.dAlpha(:,Ng+1) = a(Ng+1) .* dAlpha;
                deriv.dBeta(:)       = a(Ng+1) .* dBeta;
                deriv.dTau(:,Ng+1)   = a(Ng+1) .* dTau;

            end
        otherwise
            error('Unkown sigmoid type')
    end
end

end

