function kG = iGRAPPA(k,opt)
%   Applies the inverse grappa kernel to the undersampled k-space/image
%   
%   Author: RML
%   Date: 02/2011
%   
%   Usage: kG = iGRAPPA(k,opt)
%   
%   Input:
%   k: Undersampled k-space/image of size RO x PE x NS x NT x NR
%   opt: Specifies transform options for recon. Use SPSENSE_optset.m
%   	Must include fields:
%           opt.GKERN: GRAPPA/SPIRiT kernel of size:
%               RO x PE x NS x NT x NR x NR
%   
%   Output:
%   kG: difference signal

%   Call mex function
kG = iGRAPPA_MEX(k,opt.GKERN);

end
