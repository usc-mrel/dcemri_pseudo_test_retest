function WC = fwtN(img,wname,D)
%   Multi-level, N-dimensional forward wavelet transform
%
%   Usage:
%     WC = fwtN(img,wname)
%     WC = fwtN(img,wname,D)
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
[LPD,HPD] = get_wavelet_dec_filter(wname);

%   Convert to single precision, if necessary
if isa(img,'single');
    LPD = single(LPD);
    HPD = single(HPD);
end

%   If necessary, setup D
if nargin==2
    Nd = ndims(img);
    if size(img,Nd)==1
        Nd = Nd-1;
    end
    D = 1:Nd;
end

%   Perform transform(s)
if isnumeric(D)
    if length(D)==1 % include this case explicitly for improved 1-D speed
        [WC.A,WC.D] = fwt1dMEX(img,LPD,HPD,D);
        WC.odd_dims = isodd(size(img));
    else
        [WC] = fwtn_sub(img,LPD,HPD,D);                     
    end
    
else
    %   Preallocate the wavelet coefficient structure
    ND = length(D);
    WC = cell(1,ND);
    
    %   Loop through requested transform levels
    for ii = 1:ND        
        
        %   Perform fwd transform
        [WC{ii}] = fwtn_sub(img,LPD,HPD,D{ii});
        
        %   If necessary,
        %   prepare for next loop iteration, and wipe out the approx coeffs
        if ii<ND            
            afield = repmat('A',[1 length(D{ii})]);
            img = WC{ii}.(afield);
            WC{ii}.(afield) = [];
        end
    end
    
end

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


