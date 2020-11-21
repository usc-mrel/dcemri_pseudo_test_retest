function [S,L] = sMAP(G,opt)
%   Computes the sensitivity profile, given the GRAPPA/SPIRiT kernal
%   
%   Author: RML
%   Date: 03/2012
%   
%   Usage: [S,L] = sMAP(G,opt)
%   
%   Input:
%   kcal: Image space kernel of size RO x PE x PE2 x NT x NR x NR
%   opt: Specifies transform options for recon. Use SPSENSE_optset.m
%   
%   Output:
%   S: Dominant eigenvector map (ie, coil sensitivity)
%   L: Image mask


%   Get full size
ro = opt.size(1);
pe = opt.size(2);
ns = opt.size(3);
nr = opt.size(5);

%   Get kernal size
[Gro, Gpe, Gns, ~, ~, ~] = size(G);

%   Create output
Sl = zeros([Gro,Gpe,Gns,1,nr],opt.class);
Ll = zeros([Gro,Gpe,Gns],opt.class);

%   Loop through spatial locations and compute eigenvectors
if opt.verbose
    fprintf('Computing eigenvalues and eigenvectors\n');
end
for i = 1:Gro
    Gt = permute(G(i,:,:,:,:,:),[1 2 3 4 6 5]);
    for j = 1:Gpe
    for k = 1:Gns
        [v,l] = eig(squeeze(Gt(1,j,k,1,:,:)));
        Sl(i,j,k,1,:) = v(:,1);
        Ll(i,j,k) = l(1);
    end
    end
end
clear i j k Gt

%   Make phase relative to first receiver coil (otherwise it's pretty crazy)
Sl = Sl .* repmat(exp(-sqrt(-1)*angle(Sl(:,:,:,:,1))),[1 1 1 1 nr]);

%   Interpolate to full size
if (Gro ~= ro) || Gpe ~= pe || Gns ~= ns
    S = zeros([ro pe ns 1 nr],opt.class);
    L = zeros([ro pe ns],opt.class);
    F = fermi(ones([Gro,Gpe,Gns],opt.class),0.9,1/50);
    for i = 1:nr
        S(ceil((ro-Gro)/2)+1:ceil((ro+Gro)/2),...
          ceil((pe-Gpe)/2)+1:ceil((pe+Gpe)/2),...
          ceil((ns-Gns)/2)+1:ceil((ns+Gns)/2),...
          1,i) = F.*fFastFT(Sl(:,:,:,1,i),[1 2 3],1);
        
        S(:,:,:,1,i) = iFastFT(S(:,:,:,1,i),[1 2 3],1);
    end
    
    L(ceil((ro-Gro)/2)+1:ceil((ro+Gro)/2),...
      ceil((pe-Gpe)/2)+1:ceil((pe+Gpe)/2),...
      ceil((ns-Gns)/2)+1:ceil((ns+Gns)/2)) = fFastFT(Ll,[1 2 3],1);
    
    L = iFastFT(L,[1 2 3],1);
    L = (ro*pe*ns)/(Gro*Gpe*Gns) * L;
    
else
    S = Sl;
    L = Ll;
end

%   Perform masking
% L = abs(L);
% L(L > 1.25) = 0;
% L(L < 0.7) = 0;
% L(L > 0) = 1;
