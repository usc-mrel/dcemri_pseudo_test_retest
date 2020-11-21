function imgSM = svd_smooth_patch(img,NC,debug)

if nargin < 3 || isempty(debug)
    debug = 0;
end

%   Define block size
bl = 8;

%   Image size
[np, nv, ns, nt] = size(img);
if NC > nt
    NC=nt;
    imgSM = img;
    warning('NC >= nt; no smoothing has been performed');
    return
end

%   Pad image with zeros, if needed
if rem(np,bl) ~= 0
    img = cat(1,img,zeros([bl-rem(np,bl),nv,ns,nt]));
    NP = size(img,1);
else
    NP = np;
end
if rem(nv,bl) ~= 0
    img = cat(2,img,zeros([NP,bl-rem(nv,bl),ns,nt]));
    NV = size(img,2);
else
    NV = nv;
end
if ns > bl && rem(ns,bl) ~= 0
    img = cat(3,img,zeros([NP, NV, bl-rem(ns,bl), nt]));
    NS = size(img,3);
else
    NS = ns;
end

%   Number of blocks
nbl1 = floor(NP/bl);
nbl2 = floor(NV/bl);
nbl3 = max(floor(NS/bl),1);

%   Number of shifts
if nbl3 > 1
    shifts = [0    0    0;
              bl/2 bl/2 bl/2];
else
    shifts = [0    0    0;
              bl/2 bl/2 0];
end
imgSM = zeros([NP NV NS nt]);
for ish = 1:size(shifts,1)
    imgIN = circshift(img,shifts(ish,:));
    imgOUT = zeros([NP NV NS nt]);
    for bl1 = 1:nbl1
    for bl2 = 1:nbl2
    for bl3 = 1:nbl3
        
        ind1 = (1:bl) + (bl1-1)*bl;
        ind2 = (1:bl) + (bl2-1)*bl;
        if ns > bl
            ind3 = (1:bl) + (bl3-1)*bl;
        else
            ind3 = 1:ns;
        end
        
        imgbl = imgIN(ind1,ind2,ind3,:);
        [n1,n2,n3,n4] = size(imgbl);
        imgbl = reshape(imgbl,[n1*n2*n3 n4]);
        
        [~,S,V] = svd(imgbl,0);
        BF = V(:,1:NC);
        if debug
            figure(1);
            subplot(2,1,1);
            semilogy(diag(S),'x-');
            subplot(2,1,2);
            plot(abs(BF));
            legend(num2str((1:NC)'));
            drawnow;
        end
        clear S
        
        wt = pinv(BF)*imgbl';
        
        imgbl = (BF*wt)';
        imgOUT(ind1,ind2,ind3,:) = reshape(imgbl,[n1 n2 n3 n4]);
        
    end
    end
    end
    imgSM = imgSM + circshift(imgOUT,-shifts(ish,:));
end

imgSM = imgSM/size(shifts,1);

%   Remove zeros padding
imgSM = imgSM(1:np,1:nv,1:ns,:);

