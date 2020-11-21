function img2 = imresizeFT(img,s,fc)
%   Interpolates images in 3D

%   Check inputs
if nargin < 3
    fc = Inf;
end
if nargin < 2
    s = 2;
end
if nargin < 1
    error('Requires at least one input');
end

%   Get image size
[np, nv, ns] = size(img);
if ns ~= size(img,3);
    error('Image must be of size NP x NV x NS');
end

%   Check desired interpolation
if numel(s) == 1
    s = round(s.*[np nv ns]);
elseif numel(s) ~= 3
    error('Wrong interpolation size');
end

%   Create filter
if isstruct(fc)
    F = GEfermi(fc,[np nv]);
    F = repmat(F,[1 1 ns]);
elseif ~isinf(fc) && fc > 0
    F = fermi(ones([np nv ns],class(img)),fc,1/25);
else
    F = 1;
end
if numel(F) ~= numel(img) || any(size(F) ~= size(img))
    F = 1;
end

%   Create interpolation points
img2 = zeros([s(1) s(2) s(3)],class(img));
img2(ceil((s(1)-np)/2)+1:ceil((s(1)+np)/2),...
     ceil((s(2)-nv)/2)+1:ceil((s(2)+nv)/2),...
     ceil((s(3)-ns)/2)+1:ceil((s(3)+ns)/2)) = F.*fFastFT(img,[1 2 3],1);
img2 = iFastFT(img2,[1 2 3],1);

%   Scale img2 due to FT
img2 = (prod(s)/(np*nv*ns)) * img2;

end

function Mh = fermi(M,E,T)
%   Applies a Fermi filter to the data
%
%   Author: Marc Lebel 02/07
%   Usage:
%   Mh = fermi(M,E,T)
%   To output the filter, input a unit matrix
%
%   input:
%   M is a 1 or 2 dimensional matrix
%   E is is the "energy" (frequency cutoff) between 0 and 1 (default: 0.9)
%   T is the "temperture" (default: 1/100)
%
%   OUTPUT:
%   Mh is the filtered data of the same size

%   Check inputs
if nargin < 1
    error('fermi: function requires at least one input')
end
if nargin < 2
    E = 0.9;
end
if nargin < 3
    T = 1/25;
end

%   Size of input
sM = size(M);

if ndims(sM) > 3
    error('fermi: function can only generate 1, 2, or 3 dimentional filters')
end

%   Generate 1D filter
if sM(1) == 1 || sM(2) == 1
    x = -1:2/(length(M)-1):1;
    H = 1./(1+exp((abs(x)-E)/T));
    if size(M,1) > 1
        Mh = M.*H';
    else
        Mh = M.*H;
    end
elseif ndims(M) == 2
    %   Generate 2D filter
    x1 = -1:2/(sM(1)-1):1;
    x2 = -1:2/(sM(2)-1):1;

    H1 = 1./(1+exp((abs(x1)-E)/T))';
    H1 = repmat(H1,1,sM(2));

    H2 = 1./(1+exp((abs(x2)-E)/T));
    H2 = repmat(H2,sM(1),1);

    H = H1.* H2;
    Mh = M.*H;
elseif ndims(M) == 3
    %   Generate 3D filter (untested so needs debug)
    x1 = -1:2/(sM(1)-1):1;
    x2 = -1:2/(sM(2)-1):1;
    x3 = -1:2/(sM(3)-1):1;

    H1 = 1./(1+exp((abs(x1)-E)/T))';
    H1 = repmat(H1,[1,sM(2),sM(3)]);

    H2 = 1./(1+exp((abs(x2)-E)/T));
    H2 = repmat(H2,[sM(1),1,sM(3)]);

    H3 = 1./(1+exp((abs(x3)-E)/T));
    H3 = repmat(H3,[sM(1),1,sM(2)]);
    H3 = permute(H3,[1,3,2]);
    
    H = H1.*H2.*H3;
    Mh = M.*H;
end
end