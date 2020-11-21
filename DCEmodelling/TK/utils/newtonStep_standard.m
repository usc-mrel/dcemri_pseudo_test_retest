function [ v ] = newtonStep_standard( H, grad )
%function [ v ] = newtonStep_standard( H, grad )
%   solving
%       \theta = argmin |C(\theta) - d|^2_2
%
%   with Gauss Newton:
%       x_{k+1} = x_k - \mu Hf^-1 \gradf
%   needs to be computed.
%   This routine computes Hf^-1 \gradf for the standard model.
%   
%   Input:
%       (N: number of pixels)
%       H       Hessian matrix [3, 3, N]
%                   [ dvp^2     dKtdvp    dKepdvp ]
%                   [ dKtdvp    dKt^2     dKepdKt ]
%                   [ dKepdvp   dKepdKt   dKep^2  ]
%           
%       grad    gradient vector [3*N,1]
%                   [ dvp; dKt; dKep ]
%   Output:
%       v       new search direction
%   Yannick 2017


N = size(H,3);

dvp = grad(1:N);
dKt = grad(N+1:2*N);
dKep = grad(2*N+1:end);

d = zeros(3,N);

parfor k=1:N
   A = H(:,:,k);
   z = [dvp(k); dKt(k); dKep(k)];
   
   d(:,k) = solve_chol(A, z);
end

v = zeros(3*N,1);
v(1:N)       = d(1,:);
v(N+1:2*N)   = d(2,:);
v(2*N+1:end) = d(3,:);

end

function d = solve_chol(A, z)

   [R,p] = chol(A,'upper');   
   if p == 0
       opts.TRANSA = true;
       opts.UT     = true;
       z = linsolve(R,z,opts);
       
       opts.TRANSA = false;
       opts.UT     = true;
       d = linsolve(R,z,opts);
   elseif p == 3
       q = p-1;
       R = R(1:q,1:q);
       z = z(1:q);
       
       opts.TRANSA = true;
       opts.UT     = true;
       z = linsolve(R,z,opts);
       
       opts.TRANSA = false;
       opts.UT     = true;
       z = linsolve(R,z,opts);
       z(3) = 0;
       d = z;
   else
%        warning('Something went wrong while in voxel %i inverting the Hessian, p=%i',k,p)
        z(3) = 0;
        d = z;
        % this still needs better treatment that allows to fix the issue in
        % the next iteration of grad descent
   end
end

