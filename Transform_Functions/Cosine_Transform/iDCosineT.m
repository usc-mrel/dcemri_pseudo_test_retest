function img = iDCosineT(imgCT,opt)
%   Forward discrete cosine transform
%   
%   Author:
%   R Marc Lebel
%   10/2014
%   
%   Usage:
%   imgCT = fDCT(img,opt)

%   Get image size
[np, nv, ns, nt, nr] = size(imgCT);

if isempty(opt.DCTsz) || opt.DCTsz == 1
    img = zeros(size(imgCT),opt.class);
    for t = 1:nt
    for r = 1:nr
    for s = 1:ns
        img(:,:,s,t,r) = idct2(imgCT(:,:,s,t,r));
    end
    end
    end
    
else
    
    T = dctmtx(opt.DCTsz)';
    idct = @(block_struct) T * block_struct.data * T';
    
    img = zeros(size(imgCT),opt.class);
    for t = 1:nt
    for r = 1:nr
    for s = 1:ns
        img(:,:,s,t,r) = blockproc(imgCT(:,:,s,t,r),[opt.DCTsz opt.DCTsz],idct,'PadPartialBlocks',1);
    end
    end
    end
    
    img = img(1:opt.size(1),1:opt.size(2),:,:,:);
end

end
