function img = iCD(WCD,opt)
%   Applies the inverse reference image subtraction
%   
%   Author: RML
%   Date: 05/2012
%   
%   Usage: img = iREF(WCD,opt)
%   
%   Input:
%   WCD: Wavelet coefficient array
%   opt: Specifies transform options for recon. Use SPSENSE_optset.m
%   	Must include fields:
%           opt.IREF: reference image of size:
%               RO x PE x NS x 1 (x NR)
%           opt.class: data type
%   
%   Output:
%   img: difference image (same size)

% img = iwtN(WCD,opt.wname,opt.worder);
img = iVS(WCD,opt);
