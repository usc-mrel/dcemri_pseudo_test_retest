function WC = fcwtN(img,~,D)
%   Multi-level, N-dimensional forward wavelet transform
%
%   Usage:
%     WC = fcwtN(img,wname)
%     WC = fcwtN(img,wname,D)
%
%   Input:
%     img:        N-D array, complex or real, single or double precision
%     wname:      char array specifying type of wavelet (see get_wavelet_dec_filter)
%     D:          either an array specifying transform directions (for single level),
%                    eg: [1 2 4] for single-level 3-D transform along the dimensions 1, 2, and 4
%                 or a cell array specifying the transform directions at each level
%                    eg: {[1],[1 3],[1 2 3 4]} for 3-level transform 
%                       level 1: 1-D transform of img along dimension 1
%                       level 2: 2-D transform of level 1 result along dimensions 1 and 3
%                       level 3: 4-D transform of level 2 result along dimensions 1-4
%                 currently supports up to 4-D transforms (img can be > 4-D, though)
%                 if D is not specified, will default to single-level N-D transform
%               
%   Output:
%     WC:         structure containing wavelet coefficients (if D is an array), or
%                    cell array of wavelet coefficients (if D is a cell array)
%                 see fwtn_sub (below) for details about structure fields
%
%   Authors:
%     Travis Smith
%     R Marc Lebel
%     July, 2011
%


%   Get wavelet filters
%   Must have even length for fwt1dMEX
[Faf, ~] = FSfarras;
[af,  ~] = dualfilt1;

%   Convert to single precision, if necessary
if isa(img,'single');
    Faf{1} = single(Faf{1});
    Faf{2} = single(Faf{2});
    af{1} = single(af{1});
    af{2} = single(af{2});
end

%   If necessary, setup D
if nargin==1
    Nd = ndims(img);
    if size(img,Nd)==1
        Nd = Nd-1;
    end
    D = 1:Nd;
end

%   Preallocate the wavelet coefficient structure
ND = length(D);
WC = cell(2,ND);

%   Loop through real and imaginary tree
for ri = 1:2
    
    %   Restore image
    imgtmp = img;
    
    %   Get wavelet
    LPD(:,1) = Faf{ri}(:,1);
    HPD(:,1) = Faf{ri}(:,2);
    if length(D) > 1
        for i = 2:length(D)
            LPD(:,i) = af{ri}(:,1);
            HPD(:,i) = af{ri}(:,2);
        end
    end
    
    %   Perform transform(s)
    if isnumeric(D)
        if length(D)==1 % include this case explicitly for improved 1-D speed
            [WC{ri}.A,WC{ri}.D] = fwt1dMEX(imgtmp,LPD(:,1),HPD(:,1),D);
            WC{ri}.odd_dims = isodd(size(imgtmp));
        else
            WC{ri} = fwtn_sub(imgtmp,LPD(:,1),HPD(:,1),D);
        end
        
    else
        
        %   Loop through requested transform levels
        for ii = 1:ND
            
            %   Perform fwd transform
            [WC{ri,ii}] = fwtn_sub(imgtmp,LPD(:,ii),HPD(:,ii),D{ii});
            
            %   If necessary,
            %   prepare for next loop iteration, and wipe out the approx coeffs
            if ii<ND
                afield = repmat('A',[1 length(D{ii})]);
                imgtmp = WC{ri,ii}.(afield);
                WC{ri,ii}.(afield) = [];
            end
        end
        
    end
end

%   Combine complex coefficients
WC = CELL_cmplx_combine(WC);

end


function [WC] = fwtn_sub(img,LPD,HPD,dims)
%
%   Linear convolution followed by 2x downsampling (y = x(2:2:end))
%
%   WC: wavelet coefficient structure with fields:
%
%       odd_dims: array of length ndims(img) indicating if number of elements
%                 in each dimension of img is odd 
%                 (needed to correct ambiguity after decimation)
%
%   1-D
%       A: Approximation coefficients
%       D: Detail coefficients
%
%   2-D
%       AA = approx. vert./approx. horiz.
%       AD = approx. vert./detail  horiz.
%       DA = detail  vert./approx. horiz.
%       DD = detail  vert./detail  horiz.
%
%   3-D
%       AAA = approx. vert. / approx. horiz. / approx. through
%       AAD = approx. vert. / approx. horiz. / detail  through
%       ADA = approx. vert. / detail  horiz. / approx. through
%       ADD = approx. vert. / detail  horiz. / detail  through
%       DAA = detail  vert. / approx. horiz. / approx  through
%       DAD = detail  vert. / approx. horiz. / detail  through
%       DDA = detail  vert. / detail  horiz. / approx. through
%       DDD = detail  vert. / detail  horiz. / detail  through
%
%   4-D
%       and so on...

Nd = length(dims);

if Nd==1
    %   Take 1D transform along first direction
    [WC.A,WC.D] = fwt1dMEX(img,LPD,HPD,dims(1));
    
elseif Nd==2
    %   Take 1D transform along first direction
    [A,D] = fwt1dMEX(img,LPD,HPD,dims(1));
    
    %   Take 1D transform along second direction
    [WC.AA,WC.AD] = fwt1dMEX(A,LPD,HPD,dims(2));
    [WC.DA,WC.DD] = fwt1dMEX(D,LPD,HPD,dims(2));
    
elseif Nd==3
    %   Take 1D transform along first direction
    [A,D] = fwt1dMEX(img,LPD,HPD,dims(1));
    
    %   Take 1D transform along second direction
    [AA,AD] = fwt1dMEX(A,LPD,HPD,dims(2));
    [DA,DD] = fwt1dMEX(D,LPD,HPD,dims(2));
    
    %   Take 3D transfrom along third direction
    [WC.AAA,WC.AAD] = fwt1dMEX(AA,LPD,HPD,dims(3));
    [WC.ADA,WC.ADD] = fwt1dMEX(AD,LPD,HPD,dims(3));
    [WC.DAA,WC.DAD] = fwt1dMEX(DA,LPD,HPD,dims(3));
    [WC.DDA,WC.DDD] = fwt1dMEX(DD,LPD,HPD,dims(3));
    
elseif Nd==4
    %   Take 1D transform along first direction
    [A,D] = fwt1dMEX(img,LPD,HPD,dims(1));
    
    %   Take 1D transform along second direction
    [AA,AD] = fwt1dMEX(A,LPD,HPD,dims(2));
    [DA,DD] = fwt1dMEX(D,LPD,HPD,dims(2));
    
    %   Take 1D transfrom along third direction
    [AAA,AAD] = fwt1dMEX(AA,LPD,HPD,dims(3));
    [ADA,ADD] = fwt1dMEX(AD,LPD,HPD,dims(3));
    [DAA,DAD] = fwt1dMEX(DA,LPD,HPD,dims(3));
    [DDA,DDD] = fwt1dMEX(DD,LPD,HPD,dims(3));
    
    %   Take 1D transfrom along third direction
    [WC.AAAA,WC.AAAD] = fwt1dMEX(AAA,LPD,HPD,dims(4));
    [WC.AADA,WC.AADD] = fwt1dMEX(AAD,LPD,HPD,dims(4));
    [WC.ADAA,WC.ADAD] = fwt1dMEX(ADA,LPD,HPD,dims(4));
    [WC.ADDA,WC.ADDD] = fwt1dMEX(ADD,LPD,HPD,dims(4));
    [WC.DAAA,WC.DAAD] = fwt1dMEX(DAA,LPD,HPD,dims(4));
    [WC.DADA,WC.DADD] = fwt1dMEX(DAD,LPD,HPD,dims(4));
    [WC.DDAA,WC.DDAD] = fwt1dMEX(DDA,LPD,HPD,dims(4));
    [WC.DDDA,WC.DDDD] = fwt1dMEX(DDD,LPD,HPD,dims(4));
    
elseif Nd==0    
    WC = [];
    
else
    error('%d-D transform not supported at this time',Nd);
    
end

%   Keep track of which dimensions of the input img are odd in length
%   note that this field needs to be defined last, so that the AAA... field is
%     first in the list of structure fields (for iwtN)
WC.odd_dims = isodd(size(img));

end

function [y] = isodd(x)

y = x-2*floor(x*0.5);

end



function [af, sf] = FSfarras

% Farras filters organized for the dual-tree
% complex DWT.
%
% USAGE:
%    [af, sf] = FSfarras
% OUTPUT:
%    af{i}, i = 1,2 - analysis filters for tree i
%    sf{i}, i = 1,2 - synthesis filters for tree i
% See farras, dualtree, dualfilt1.
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

af{1} = [
                  0                  0
  -0.08838834764832  -0.01122679215254
   0.08838834764832   0.01122679215254
   0.69587998903400   0.08838834764832
   0.69587998903400   0.08838834764832
   0.08838834764832  -0.69587998903400
  -0.08838834764832   0.69587998903400
   0.01122679215254  -0.08838834764832
   0.01122679215254  -0.08838834764832
                  0                  0
 ];
   
sf{1} = af{1}(end:-1:1, :);

af{2} = [
   0.01122679215254                  0
   0.01122679215254                  0
  -0.08838834764832  -0.08838834764832
   0.08838834764832  -0.08838834764832
   0.69587998903400   0.69587998903400
   0.69587998903400  -0.69587998903400
   0.08838834764832   0.08838834764832
  -0.08838834764832   0.08838834764832
                  0   0.01122679215254
                  0  -0.01122679215254
];

sf{2} = af{2}(end:-1:1, :);
end


function [af, sf] = dualfilt1

% Kingsbury Q-filters for the dual-tree complex DWT
%
% USAGE:
%    [af, sf] = dualfilt1
% OUTPUT:
%    af{i}, i = 1,2 - analysis filters for tree i
%    sf{i}, i = 1,2 - synthesis filters for tree i
%    note: af{2} is the reverse of af{1}
% REFERENCE:
%    N. G. Kingsbury,  "A dual-tree complex wavelet
%    transform with improved orthogonality and symmetry
%    properties", Proceedings of the IEEE Int. Conf. on
%    Image Proc. (ICIP), 2000
% See dualtree
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

% These cofficients are rounded to 8 decimal places.

af{1} = [
   0.03516384000000                  0
                  0                  0
  -0.08832942000000  -0.11430184000000
   0.23389032000000                  0
   0.76027237000000   0.58751830000000
   0.58751830000000  -0.76027237000000
                  0   0.23389032000000
  -0.11430184000000   0.08832942000000
                  0                  0
                  0  -0.03516384000000
 ];
 
af{2} = [
                  0  -0.03516384000000
                  0                  0
  -0.11430184000000   0.08832942000000
                  0   0.23389032000000
   0.58751830000000  -0.76027237000000
   0.76027237000000   0.58751830000000
   0.23389032000000                  0
  -0.08832942000000  -0.11430184000000
                  0                  0
   0.03516384000000                  0
];
 
sf{1} = af{1}(end:-1:1, :);
 
sf{2} = af{2}(end:-1:1, :);
end


function WC1 = CELL_cmplx_combine(WC)

%   Determine transform order and extract the number of frames at each one
[~,N] = size(WC);

%   Loop through transform order and the various coefficients
WC1 = WC;
cnum = sqrt(-1);
for i = 1:N
    WCtr = WC{1,i};
    WCti = WC{2,i};
    fnames = fieldnames(WCtr);
    nf = size(fnames);
    for j = 1:nf
        fname = char(fnames(j));
        if ~strcmp(fname,'odd_dims')% && ~isempty(WCt1.(fname)) && ~isempty(WCt2.(fname))
            
            if ~isempty(WCtr.(fname))
                coefs_r = 0.5*WCtr.(fname);
                coefs_i = (0.5*cnum)*WCti.(fname);
                WCtr.(fname) = coefs_r + coefs_i;
                WCti.(fname) = coefs_r - coefs_i;
            end
        end
    end
    WC1{1,i} = WCtr;
    WC1{2,i} = WCti;
end

end