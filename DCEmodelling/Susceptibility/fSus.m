function deltaB = fSus(C, opt)
%function deltaB = fSus(C, opt)
%   forward susceptibility operator
%
%   Inputs
%       C       concentration of contrast agent [nx, ny, nz, nt]
%       opt     option container
%           .SIdirection    dimension of S/I direction
%           .     voxel size
%           .TK.chi         susceptibility of contrast agent
%           .TK.Hct         Hematocrit value
%
%   Output
%       deltaB  difference in B field
%
%   Example usage
%       C = randn(64, 64, 64, 50);
%       opt.SIdirection = 1;            % B0/z direction on first dimension
%       opt.voxel_size = [2, 1, 1];     % coarse resolution in z direction
%       opt.TK.chi     = 320e-9;        % susceptibility of contrast agent [2]
%
%       dB = fSus(C, opt);
%
%   References:
%        1. Bagher-Ebadian H, Jain R, Nejad-Davarani SP, et al. Model selection for DCE-T1 studies in glioblastoma. Magn. Reson. Med. 2012;68:241?251 doi: 10.1002/mrm.23211.
%        2. J Korporaal, et al. "Phase-based arterial input function measurements in the femoral arteries for quantification of dynamic contrast-enhanced (DCE) MRI and comparison with DCE-CT," MRM, 2011, 1267-1274.
%
%
%   Yannick Bliesener 2020
%

perm = circshift(1:3, [0, 3-opt.SIdirection]);
perm = [perm 4];

C = permute(C, perm);
voxelsize = opt.voxel_size(perm(1:3));

[Nx,Ny,Nz,Nt] = size(C);

% Careful with this, if TK parameter for Chi and Hct ain't the same in vessel 
% and tissue this is wrong
C = C * opt.TK.chi;

Dkernel = dipole_kernel([Nx Ny Nz], voxelsize, 3, 'kspace', 'dipole');

deltaB = ifftnd(repmat(Dkernel, [1 1 1 Nt]) .* fftnd(C, [1 2 3], 0), [1 2 3], 0);
deltaB = real(deltaB);

deltaB = ipermute(deltaB, perm);

end

