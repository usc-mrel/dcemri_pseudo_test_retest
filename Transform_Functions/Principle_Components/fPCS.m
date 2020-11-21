function PCwt = fPCS(img,opt)
%   Computes weights for principle components

%   Get image size then reshape
[np, nv, ns, nt] = size(img);
img2 = reshape(img,[np*nv*ns nt])';

%   Get basis functions
BF = opt.BF;
[ntBF,nBF] = size(BF);

if nBF > ntBF || ntBF ~= nt
    error('fPCS: wrong sizes');
end

%   Get weights
PCwt = pinv(BF)*img2;

end

