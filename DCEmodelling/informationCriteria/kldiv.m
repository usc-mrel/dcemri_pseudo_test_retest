function [ kl ] = kldiv( P, Q )
%[ kl ] = kldiv( P, Q )
%   Kullback Leibler Divergence
%
%   Yannick Bliesener 2017
%

if any( (abs(Q) < eps) & (abs(P) > eps) )
   warning( 'KL divergence is defined only if for all i: Q(i)=0 => P(i)=0' ) 
end
kl = P .* log( P ./ (Q+eps) );
kl( abs(P) < eps ) = 0;
kl = sum( kl, 1 );

end

