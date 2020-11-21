function opt = SPSENSE_optset
%   Set options for conjugate gradient SPSENSE reconstruction
%   
%   Author: RML
%   Date: 09/2011
%   
%   Usage: opt = SPSENSE_optset
%   
%   Input:
%   (none)
%   
%   Output:
%   opt: options structure

%   Define general options
opt.MaxIter = 300;              %   Maximum conjugate-gradient iterations
opt.MinIter = 1;                %   Minimum conjugate-gradient iterations
opt.class = 'single';           %   Data class for computation (double or single)
opt.sqrt_smooth = 1e-5;         %   Smoothing factor for abs() approximation
opt.Dobj = 1e-3;                %   Stopping threshold (% chanage in objective)
opt.Dobjn = 6;                 %   Median window for changes below threshold
opt.plot = 0;                   %   Flag for plotting
opt.plotf = @(img,opt)recon_plot(img,opt);  %   Plotting function
opt.verbose = 2;                %   Flag for text output (0: none; 1: some; 2: lots)
opt.update = 5;                 %   Periodic update interval
opt.compression_levels = 0;     %   Auto adjust constraints
opt.target_CL = 90;             %   Target compression level (% of coeff. below 0.05)
opt.size = [];                  %   Image size (placeholder)
opt.U = [];                     %   Sampled indices in k-space

%   Define phase correction options
opt.phase_corr = [];           %   Phase correction dimensions
opt.Ph_cf = [0.25 0.05 0.05];       %   Cutoff frequencies for phase estimation
opt.ePH = @(img,opt)estPH(img,opt); %   Phase estimation function
opt.fPH = @(img,opt)fPH(img,opt);   %   Forward phase correction
opt.iPH = @(img,opt)iPH(img,opt);   %   Inverse phase correction

%   Define FT options
opt.fFT = @(f,opt)fFastFT(f,opt);   %   Define the forward Fourier operator
opt.iFT = @(f,opt)iFastFT(f,opt);   %   Define the inverse Fourier operator
opt.fFTU= @(f,opt)fFFTU(f,opt);     %   Define the forward Fourier operator (undersampled)
opt.iFTU= @(f,opt)iFFTU(f,opt);     %   Define the inverse Fourier operator (undersampled)
opt.FTshift = 0;                    %   FT shift
opt.FTdim = [2 3];                  %   FT dimension
opt.kW = [];                        %   Optional weighting vector for k-space

%   Define GRAPPA/SPIRiT options
opt.lambda(1) = 0;                  %   Lagrange multiplier for GRAPPA penalty
opt.compression_levels(1) = 0;      %   Compression level
% opt.fPI = @(img,opt)fGRAPPA(img,opt);   %   Forward Parallel imaging function
% opt.iPI = @(img,opt)iGRAPPA(img,opt);   %   Inverse Parallel imaging function
opt.fPI = @(img,opt)fSENSE(img,opt);    %   Forward SENSE transform
opt.iPI = @(img,opt)iSENSE(img,opt);    %   Inverse SENSE transform
opt.kernreg = 'auto';               %   Regularizing factor while solving grappa kernel
opt.kernel = [3 7 7];               %   Kernel size
opt.gKERN = @(kcal,opt)gKERN3d(kcal,opt);   %   GRAPPA kernel computation

%   Define wavelet options
opt.lambda(2) = 1e-3;           %   Lagrange multiplier for wavelet penalty
opt.compression_levels(2) = 50; %   Compression level
opt.wname = 'db3';              %   Wavelet name (not used with complex wavelet transform)
opt.fSP1 = @(img,wname,worder)fcwtN(img,wname,worder);  %   Wavelet transform
opt.iSP1 = @(WC,wname,worder)icwtN(WC,wname,worder);    %   Inverse wavelet transform
opt.worder = {'auto'};          %   Placeholder (defined below)

%   Define total variation options
opt.lambda(3) = 1e-3;           %   Lagrange multiplier for finite difference penalty
opt.compression_levels(3) = 25; %   Compression level
opt.fSP2 = @(x,order)fFD(x,order);      %   Forward F.D. function
opt.iSP2 = @(x,order)iFD(x,order);      %   Inverse F.D. function
opt.FDorder1 = {1, 2};          %   Placeholder (defined below)

%   Define total variation options
opt.lambda(4) = 0;              %   Lagrange multiplier for spatial filter
opt.compression_levels(4) = 25; %   Compression level
opt.fSP3 = @(x,opt)fImFilt(x,opt);      %   Forward transform
opt.iSP3 = @(x,opt)iImFilt(x,opt);      %   Inverse transform
opt.SF_kern = [];               %   Placeholder (defined below)

%   Define view sharing options
opt.lambda(5) = 5e-3;           %   Lagrange multiplier for view-sharing (for dynamic scans)
opt.compression_levels(5) = 80; %   Compression level
opt.vs_kern = 1/sqrt(2)*[-1 1 0];%   View-sharing kernel
opt.fVS = @(x,opt)fVS(x,opt);   %   Forward view sharing function
opt.iVS = @(x,opt)iVS(x,opt);   %   Inverse view sharing funciton

%   Define reference image options
opt.lambda(6) = 1e-3;           %   Lagrange multiplier for difference relative to reference image
opt.compression_levels(6) = 50; %   Compression level
opt.fSP4 = @(x,opt)fREF(x,opt); %   Forward transform function
opt.iSP4 = @(x,opt)iREF(x,opt); %   Inverse transform function
opt.IREF = [];                  %   Placeholder for reference image
opt.Npc = 2;                    %   Number of principle components for smoothing
opt.Npc_Nup = 50;               %   Number of updates to perform
opt.BF = [];                    %   Optional basis functions for PCA
opt.Npc_local = 1;              %   Perform global PCA or local patch based

%   Define spatial support
opt.lambda(7) = 1e-3;           %   Lagrange multiplier for difference relative to reference image
opt.compression_levels(7) = 10; %   Compression level
opt.fCD = @(x,opt)fCD(x,opt);   %   Forward transform function
opt.iCD = @(x,opt)iCD(x,opt);   %   Inverse transform function

%   Define l1-weighting options
opt.l1weight = 0;               %   Exponent to: 1/(X^l1weight + l1reg)
opt.l1reg = 1e-2;               %   Smoothing factor to prevent Inf weight
opt.fWT = @(WCwt,WC)fWT(WCwt,WC);   %   Weighting functions
