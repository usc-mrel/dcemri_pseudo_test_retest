function imgCT = fDCosineT(img,opt)
%   Forward discrete cosine transform
%   
%   Author:
%   R Marc Lebel
%   10/2014
%   
%   Usage:
%   imgCT = fDCT(img,opt)

%   Get image size
[np, nv, ns, nt, nr] = size(img);

if isempty(opt.DCTsz) || opt.DCTsz == 1
    imgCT = zeros(size(img),opt.class);
    for t = 1:nt
    for r = 1:nr
    for s = 1:ns
        imgCT(:,:,s,t,r) = dct2(img(:,:,s,t,r));
    end
    end
    end
    
else
    T = dctmtx(opt.DCTsz);
    dct = @(block_struct) T * block_struct.data * T';
    
    sze = opt.DCTsz*ceil([np,nv]/opt.DCTsz);
    imgCT = zeros([sze ns nt nr],opt.class);
    for t = 1:nt
    for r = 1:nr
    for s = 1:ns
        imgCT(:,:,s,t,r) = blockproc(img(:,:,s,t,r),[opt.DCTsz opt.DCTsz],dct,'PadPartialBlocks',1);
    end
    end
    end
end

end
