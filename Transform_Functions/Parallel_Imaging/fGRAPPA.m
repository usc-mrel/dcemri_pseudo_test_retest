function kG = fGRAPPA(k,opt)
%   Applies the forward grappa kernel to the undersampled k-space/image
%   
%   Author: RML
%   Date: 02/2011
%   
%   Usage: kG = fGRAPPA(k,opt)
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
kG = fGRAPPA_MEX(k,opt.GKERN);

end