function img = iPCS(PCwt,opt)
%   Computes weights for principle components


%   Get basis functions
BF = opt.BF;
[ntBF,nBF] = size(BF);

if nBF > ntBF || ntBF ~= opt.size(4)
    error('iPCS: wrong sizes');
end

%   Apply weights to get image
img = BF*PCwt;
img = img';
img = reshape(img,opt.size(1:4));

end

