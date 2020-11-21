function [ Y, S ] = truncateSingValMatrix( X, k )
%function [ Y, S ] = truncateSingValMatrix( X, k )
%   computes low rank approximation of a matrix
%
% Yannick 2017

[m, n] = size(X);

[U,S,V] = svd(X, 'econ');

cutOff = min([m, n, k]);

if (size(S,2) == 1) || (size(S,1) == 1)
    S = S(1,1);
else
	S = diag( S );
end
Scut = S(1:cutOff);

U = U(:,1:cutOff);
V = V(:,1:cutOff);

Y = U*diag(Scut)*V';

end

