function img = fDESPOT1(T1,Mo,opt)
%   Converts R1 and Mo maps into images
%   
%   opt: options structure with fields
%       FA: flip angle (degrees)
%       tr: repetition time (s)
%       B1: transit field scale (fraction)

[np,nv,ns] = size(T1);
nFA = length(opt.FA);

%   B1
if ~isfield(opt,'B1') || isempty(opt.B1)
    opt.B1 = ones([np, nv, ns]);
end

%   Use SPGR equation to compute images
img = zeros([np,nv,ns,nFA]);
E1 = exp(-opt.tr./T1);
E1m1 = (1-E1);
for iFA = 1:nFA
    fa = opt.B1 .* opt.FA(iFA) * pi/180;
    img(:,:,:,iFA) = Mo .* (E1m1).*sin(fa) ./ (1 - E1.*cos(fa));
end

end
