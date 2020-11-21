function img = icwtN(WC,~,D)
%   Multi-level multi-dimensional inverse wavelet transform
%
%   Usage:
%     img = icwtN(WC,wname)
%     img = icwtN(WC,wname,D)
%
%   Input:
%     WC:         output from fwtN_
%                 complex or real, single or double precision
%     wname:      char array specifying type of wavelet (see get_wavelet_rec_filter)
%     D:          the array specifying the levels/dimension which was used to create WC (see fwtN)
%                 currently supports up to 4-D transforms (img can be > 4-D)
%                 if not specified, will default to single-level N-D transform
%
%   Output:
%     img:        N-D array containing the result
%
%   Authors:
%     Travis Smith
%     R Marc Lebel
%     July, 2011
%
%

%   Get wavelet filters
%   Must have even length for fwt1dMEX
[~, Fsf] = FSfarras;
[~,  sf] = dualfilt1;

%   Convert to single precision, if necessary
if isstruct(WC)
    fnames = fieldnames(WC);
    s = isa(WC.(fnames{1}),'single');
else
    fnames = fieldnames(WC{end});
    s = isa(WC{end}.(fnames{1}),'single');
end

if s
    Fsf{1} = single(Fsf{1});
    Fsf{2} = single(Fsf{2});
    sf{1} = single(sf{1});
    sf{2} = single(sf{2});
end


%   If necessary, setup D
if nargin==2 && isstruct(WC)
    Nd = log2(length(fieldnames(WC))-1);
    D = 1:Nd;
end

%   Undo complex dual tree coefficients
WC = CELL_cmplx_extract(WC);

%   Loop through real and imaginary tree
imgtmp = cell(1,2);
for ri = 1:2
    
    %   Get wavelet
    LPR(:,1) = Fsf{ri}(:,1);
    HPR(:,1) = Fsf{ri}(:,2);
    if length(D) > 1
        for i = 2:length(D)
            LPR(:,i) = sf{ri}(:,1);
            HPR(:,i) = sf{ri}(:,2);
        end
    end
    
    %   Perform transform(s)
    if isnumeric(D)
        if length(D)==1 % include this case explicitly for improved 1-D speed
            img = iwt1dMEX(WC{ri}.A,WC{ri}.D,LPR(:,1),HPR(:,1),D,WC{ri}.odd_dims(D));
        else
            img = iwtn_sub(WC{ri},LPR(:,1),HPR(:,1),D);
        end
        
    else
        ND = length(D);
        fnames = fieldnames(WC{ri,end});
        img = WC{ri,ND}.(fnames{1}); %  This assmumes that the approximation coefficient is first...
        
        %   Loop through requested transform levels (in reverse)
        for ii = ND:-1:1
            
            fnames = fieldnames(WC{ri,ii});
            if length(fnames)==1 % 1 field name is for input_size
                if ~isempty(D{ii})
                    error('User specified %d-D transform, but no wavelet coeffs found',length(D{ii}));
                else
                    continue;
                end
            end
            
            %   If necessary,
            %   prepare for transform
            if ii < ND
                WC{ri,ii}.(fnames{1}) = img; %  This assmumes that the approximation coefficient is first...
            end
            
            %   Perform inv transform
            img = iwtn_sub(WC{ri,ii},LPR(:,ii),HPR(:,ii),D{ii});
            
        end
        
    end
    
    imgtmp{ri} = img;
    
end

%   Combine image
img = imgtmp{1} + imgtmp{2};

end


function img = iwtn_sub(WC,LPR,HPR,dims)
%
%   2x upsampling (y(2:2:end) = x) followed by linear convolution and cropping
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
    img = iwt1dMEX(WC.A,WC.D,LPR,HPR,dims(1),WC.odd_dims(dims(1)));
    
elseif Nd==2
    %   Combine second direction coefficients
    p = WC.odd_dims(dims(2));
    A = iwt1dMEX(WC.AA,WC.AD,LPR,HPR,dims(2),p);
    D = iwt1dMEX(WC.DA,WC.DD,LPR,HPR,dims(2),p);
    
    %   Combine first direction coefficients
    img = iwt1dMEX(A,D,LPR,HPR,dims(1),WC.odd_dims(dims(1)));
    
elseif Nd==3
    %   Combine third direction coefficients
    p = WC.odd_dims(dims(3));
    AA = iwt1dMEX(WC.AAA,WC.AAD,LPR,HPR,dims(3),p);
    AD = iwt1dMEX(WC.ADA,WC.ADD,LPR,HPR,dims(3),p);
    DA = iwt1dMEX(WC.DAA,WC.DAD,LPR,HPR,dims(3),p);
    DD = iwt1dMEX(WC.DDA,WC.DDD,LPR,HPR,dims(3),p);
    
    %   Combine second direction coefficients
    p = WC.odd_dims(dims(2));
    A = iwt1dMEX(AA,AD,LPR,HPR,dims(2),p);
    D = iwt1dMEX(DA,DD,LPR,HPR,dims(2),p);
        
    %   Combine first direction coefficients
    img = iwt1dMEX(A,D,LPR,HPR,dims(1),WC.odd_dims(dims(1)));
    
elseif Nd==4
    %   Combine fourth dimension coefficients
    p = WC.odd_dims(dims(4));
    AAA = iwt1dMEX(WC.AAAA,WC.AAAD,LPR,HPR,dims(4),p);
    AAD = iwt1dMEX(WC.AADA,WC.AADD,LPR,HPR,dims(4),p);
    ADA = iwt1dMEX(WC.ADAA,WC.ADAD,LPR,HPR,dims(4),p);
    ADD = iwt1dMEX(WC.ADDA,WC.ADDD,LPR,HPR,dims(4),p);
    DAA = iwt1dMEX(WC.DAAA,WC.DAAD,LPR,HPR,dims(4),p);
    DAD = iwt1dMEX(WC.DADA,WC.DADD,LPR,HPR,dims(4),p);
    DDA = iwt1dMEX(WC.DDAA,WC.DDAD,LPR,HPR,dims(4),p);
    DDD = iwt1dMEX(WC.DDDA,WC.DDDD,LPR,HPR,dims(4),p);
    
    %   Combine third direction coefficients
    p = WC.odd_dims(dims(3));
    AA = iwt1dMEX(AAA,AAD,LPR,HPR,dims(3),p);
    AD = iwt1dMEX(ADA,ADD,LPR,HPR,dims(3),p);
    DA = iwt1dMEX(DAA,DAD,LPR,HPR,dims(3),p);
    DD = iwt1dMEX(DDA,DDD,LPR,HPR,dims(3),p);
    
    %   Combine second direction coefficients
    p = WC.odd_dims(dims(2));
    A = iwt1dMEX(AA,AD,LPR,HPR,dims(2),p);
    D = iwt1dMEX(DA,DD,LPR,HPR,dims(2),p);
    
    %   Combine first direction coefficients
    img = iwt1dMEX(A,D,LPR,HPR,dims(1),WC.odd_dims(dims(1)));
    
else
    error('%d-D transform not supported at this time',Nd);
end

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

function WC1 = CELL_cmplx_extract(WC)

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
                coefs_i = 0.5*WCti.(fname);
                WCtr.(fname) = coefs_r + coefs_i;
                WCti.(fname) = -cnum*coefs_r + cnum*coefs_i;
                WCtr.(fname)(1) = eps+eps*cnum;
                WCti.(fname)(1) = eps+eps*cnum;
            end
        end
    end
    WC1{1,i} = WCtr;
    WC1{2,i} = WCti;
end

end

