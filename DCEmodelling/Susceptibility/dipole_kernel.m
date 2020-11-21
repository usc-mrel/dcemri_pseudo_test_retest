function D = dipole_kernel( matrix_size, voxel_size, B0_dir, idomain, ishape, iaug )
%function D = dipole_kernel( matrix_size, voxel_size, B0_dir, idomain, ishape, iaug )
% Generation of the Dipole Kernel
%
%   Input
%       matrix_size     the size of the matrix
%       voxel_size      the size of a voxel
%       B0_dir          the direction of the B0 field
%       idomain         'imagespace' or 'kspace'
%       ishape          elementary shape of kernel
%       iaug            augmentation factor to augment dimension of elementary shape across k-space grid                
%           
%   Output
%       D       dipole kernel saved in Fourier space
% 
%   Refernces:
%       Fourier domain expression:
%           Salomir et al. Concepts in Magn Reson Part B 2003, 19(B):26-34
%       Image domain expression:
%           Li et al. Magn Reson Med 2004, 51(5):1077-82
%
%   Created by Tian Liu in 2008
%   Modified by Tian Liu on 2011.02.01
%   Last modified by Tian Liu on 2013.07.22
%   Modified by Yannick Bliesener 2017
%
%   TODO:
%       - maintain imagespace
%       - make code ready for other B0 directions than z

%% Parse input
lorentzCorr     = true;
    % apply Lorentz correction
softEquivalence = false;

domain      = 'kspace';
    % domain of kernel
shape       = 'dipole';
    % elementary shape of kernel
aug         = 1;
    % augmentation factor to augment the k-space grid in order to resolve
    % kernel shapes other than simple dipoles

if nargin > 3
    domain = idomain;
end
if nargin > 4
    shape = ishape;
end

if nargin > 5
    aug = iaug;
end

if B0_dir ~= 3
    error('Source code currently not ready for any other direction than z!')
end

if (B0_dir == 1)
    B0_dir = [1 0 0 ]';
elseif (B0_dir == 2)
    B0_dir = [0 1 0 ]';
elseif (B0_dir==3)
    B0_dir = [0 0 1]';
end


%% Compute dipole kernel

switch domain
    case 'kspace'
        [Y,X,Z]=meshgrid(-matrix_size(2)/2:(matrix_size(2)/2-1),...
            -matrix_size(1)/2:(matrix_size(1)/2-1),...
            -matrix_size(3)/2:(matrix_size(3)/2-1));
        
        % kspace normalized from -0.5 to 0.5
        X = X/(matrix_size(1)*voxel_size(1));
        Y = Y/(matrix_size(2)*voxel_size(2));
        Z = Z/(matrix_size(3)*voxel_size(3));
        
        % kspace normalized from -1 to 1
%         X = X/(max(abs(X(:)))*voxel_size(1));
%         Y = Y/(max(abs(Y(:)))*voxel_size(2));
%         Z = Z/(max(abs(Z(:)))*voxel_size(3));


        % Gaussian quadrature
        if softEquivalence
            gauss = 0;
            
            quadGrid    = [-sqrt(1/3) +sqrt(1/3)] / 2;
            quadWeight  = [1 1] / 2;
            
            for i=1:2
                for j=1:2
                    for k=1:2
                        tmp = quadWeight(i)*quadWeight(j)*quadWeight(k);
                        tmp = tmp .* exp(-1i*2*pi * X * voxel_size(1) * quadGrid(k));
                        tmp = tmp .* exp(-1i*2*pi * Y * voxel_size(2) * quadGrid(j));
                        tmp = tmp .* exp(-1i*2*pi * Z * voxel_size(3) * quadGrid(i));
                        gauss = gauss + tmp;
                    end
                end
            end
        else
            gauss = 1;
        end
        
        if lorentzCorr
            % Lorentz corrected:
            D = 1/3-  ( X*B0_dir(1) + Y*B0_dir(2) + Z*B0_dir(3) ).^2./(X.^2+Y.^2+Z.^2);
        else
            % without correction:
            D = 1-  ( X*B0_dir(1) + Y*B0_dir(2) + Z*B0_dir(3) ).^2./(X.^2+Y.^2+Z.^2);
        end
        
        switch shape
            case {'dipole'}
                % D = D .* 1;

                D = D .* gauss;
            case {'box'}
                V = 1;
                a = aug*voxel_size(1);
                b = aug*voxel_size(2);
                c = aug*voxel_size(3);
                
                D = D .* V .* sinc(a*X/(2*pi)) .* sinc(b*Y/(2*pi)) .* sinc(c*Z/(2*pi));
            case {'sphere'}
%                 V = prod(voxel_size);
                V = 1;
                R = aug*min(voxel_size(1:2));
%                 R = nthroot(3*V/4/pi,3);
                kR = R*sqrt( X.^2 + Y.^2 + Z.^2);

                D = D .* (3*V) .* spherical_besselj(1,kR) ./ kR;
                
            case 'cylinder'
                h = aug*voxel_size(3);
                R = aug*min(voxel_size(1:2));
%                 V = pi*h*R^2;
%                 V = prod(voxel_size);
                V = 1;
                
                D = D .* V .* sinc(h*Z/(2*pi)) .* jinc(R*sqrt( X.^2 + Y.^2)/pi); % .* exp(-1j*pi*Y);
            otherwise
                error('Unkown shape for kspace domain')
        end

        D(isnan(D)) = 0;
        D = fftshift(D);
        
    case 'imagespace'
        [Y,X,Z]=meshgrid(-matrix_size(2)/2:(matrix_size(2)/2-1),...
            -matrix_size(1)/2:(matrix_size(1)/2-1),...
            -matrix_size(3)/2:(matrix_size(3)/2-1));
        
        X = X*voxel_size(1);
        Y = Y*voxel_size(2);
        Z = Z*voxel_size(3);
        
        switch shape
            case {'sphere', 'dipole'}
                
                d = prod(voxel_size)*(3*( X*B0_dir(1) + Y*B0_dir(2) + Z*B0_dir(3)).^2 - X.^2-Y.^2-Z.^2)./(4*pi*(X.^2+Y.^2+Z.^2).^2.5);
                d(isnan(d)) = 0;
                
            case 'infcylinder'
                d = zeros( size(X) );
                d( sqrt(X.^2 + Y.^2 + Z.^2) <= eps ) = 1/3;
            otherwise
                error('Unkown shape for imagespace domain')
        end
        
        D = fftn(fftshift(d));
    otherwise
        error( 'Unkown domain!' )
end


end