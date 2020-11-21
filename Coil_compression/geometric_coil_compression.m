function [k2, nr2] = geometric_coil_compression(k, nr2, petab)
%   Geometric coil compression

[np, nvtmp, nstmp, nttmp, nr] = size(k);

%   Find good calibration points (high signal in k-space)
if nargin == 3 && ~isempty(petab)
    ncal = 8192;
    nv = max(petab(:,1)) + min(petab(:,1)) - 1;
    ns = max(petab(:,2)) + min(petab(:,2)) - 1;
    calind = find(petab(:,1)>nv/2-15 & petab(:,1)<nv/2+16 & ...
        petab(:,2)>ns/2-15 & petab(:,2)<ns/2+16,ncal);
    kcal = k(:,calind,:,:);
else
    %   Assume image with [readout x phase x phase2 x other x coils]
    kcal = squeeze(k(:,:,:,1,:));
    
end

%   Compute unaligned compression matrix
gccmtx = calcGCCMtx(kcal,1,5);

%   Determine a suitable number of coils
%   Done by looking at the k-space power after compression
if isempty(nr2) || nr2 == 0
%     kcal2 = CC(kcal,1,gccmtx);
%     l2nrm_1 = zeros(1,nr);
%     l2nrm_2 = zeros(1,nr);
%     for i = 1:nr
%         l2nrm_1(i) = norm(kcal(:,:,:,i));
%         l2nrm_2(i) = norm(kcal2(:,:,:,i));
%     end
%     
%     %nr2 = sum(l2nrm_2 > mean(l2nrm_1)/10);
%     nr2 = sum(l2nrm_2 > min(l2nrm_1)/10);

    nr2 = ceil(18.*(1-exp(-nr/30)));
end

%   Compute aligned matrix
if nr2 <= nr
    gccmtx = alignCCMtx(gccmtx(:,1:nr2,:));
end

%   Apply to all data
k2 = CC(reshape(k,[np,nvtmp*nstmp*nttmp,1,nr]),1,gccmtx);
k2 = reshape(k2,[np, nvtmp, nstmp, nttmp, nr2]);

end

