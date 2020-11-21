%%
%   This script is the central reconstruction script of the STARDCE
%   package.
%   It is supposed to be executed at the end of the DCE exam.
%
%
%
%
%
%

clear all; close all;

%% load dependencies
addpath( genpath('/server/home/dce/orchestra-sdk-1.8-1.matlab') );
addpath( genpath('/server/home/dce/STARDCE/T1estimation') );
addpath( genpath('/server/home/dce/Gpufit/Gpufit-build/Matlab') );

% to plot results
addpath( genpath('/server/home/dce/export_fig') );
addpath( genpath('/server/home/dce/tight_subplot') );

%% determine folder to reconstruct
path_to_data = getenv('STARDCE_DATA_PATH');

cd(path_to_data)

%% set global parameters

matrixSize = [256 240 120];

%% 1) compute T1 maps
GRopt = GRoptset(8, 't1');

GRopt.DO_SPSENSE = 0;    
GRopt.DO_DIRECT  = 1;   

GRopt.mtrxH = matrixSize;

GRopt.bolus_arrival = 300;
    
[T1, Mo, ~] = GRT1map(GRopt);

%% 2) reconstruct pre and post contrast (for brain mask)
GRopt = GRoptset(64, 'perm');

GRopt.mtrxH = matrixSize;

GRopt.DO_MOCCO = 0;
GRopt.DO_PRE_POST = 1;

GRopt.bolus_arrival = 300;

[~, imgPRE, ~, imgPOST, ~, ~, ~] = GRrecon(GRopt, [], 'pre_post');

%% 3) reconstruct dynamic images
GRopt = GRoptset(64, 'perm');

GRopt.mtrxH = matrixSize;

GRopt.DO_MOCCO = 0;
GRopt.DO_SPSENSE = 1;

GRopt.bolus_arrival = 300;

[~, ~, ~, ~, imgR, OUT] = GRrecon(GRopt, [], 'SPSENSE');

%% 4) estimate TK parameters

% sparse SENSE based TK estimation
[~, TK, TKESTOUT] = GRestimation(GRopt, [], 'SPSENSE');

% MOCCO based estimation
% GRopt = GRoptset(64, 'perm');
% 
% GRopt.mtrxH = matrixSize;
%     
% [imgR, TK, OUT] = GRestimation(GRopt, path, 'MOCCO');

%% 5) make output figures
close all;

saggital = 112;
coronal = 55;
axial = 178;

nrows = 5;
ncols = 3;

brainmask = ones(matrixSize);

% %%%%%%%%%%%%%% center-of-kspace
load('kcent.mat')
h1 = figure(1);
plot(tcent, kcent, 'LineWidth', 2)
xlabel('Time [s]')
ylabel('Signal intensity of k-space center')
title('Center of k-space')
set(h1, 'Color', 'w');
set(findall(h1, '-property','FontSize'),'FontSize',18)
export_fig(h1, 'kcent.pdf', '-p0.02')



%   in:  Nh      number of axes in hight (vertical direction)
%        Nw      number of axes in width (horizontaldirection)
%        gap     gaps between the axes in normalized units (0...1)
%                   or [gap_h gap_w] for different gaps in height and width 
%        marg_h  margins in height in normalized units (0...1)
%                   or [lower upper] for different lower and upper margins 
%        marg_w  margins in width in normalized units (0...1)
%                   or [left right] for different left and right margins 

h2 = figure(2);
hax = tight_subplot(nrows, ncols, [0. 0], [0.05 0.05], [0.2 0.2]);
hax = reshape(hax, [ncols, nrows])';

% %%%%%%%%%%%%%% B1+ (axial, sagittal, coronal)
load('GRrecon_T1Mo.mat', 'B1')

axes(hax(1, 1));
im = permute(B1(:,:,:) .* brainmask, [3 2 1]);
im = flip(im,1);
imshow(im(:,:,axial), [0.5 1.5])
axis square
colormap(hax(1, 1), 'jet')
ylabel('B1+', 'rotation',0,'VerticalAlignment','middle');

axes(hax(1, 2));
im = permute(B1(:,:,:) .* brainmask, [1 2 3]);
imshow(im(:,:,coronal), [0.5 1.5])
axis square
colormap(hax(1, 2), 'jet')

axes(hax(1, 3));
im = permute(B1(:,:,:) .* brainmask, [1 3 2]);
im = flip(im, 2);
imshow(im(:,:,saggital), [0.5 1.5])
axis square
colormap(hax(1, 3), 'jet')
pos = get(hax(1, 3), 'Position');
hc = colorbar;
ylabel(hc, '[unitless]')
set(hax(1, 3), 'Position', pos)

% %%%%%%%%%%%%%% M0 (axial, sagittal, coronal)
load('GRrecon_T1Mo.mat', 'Mo')

axes(hax(2, 1));
im = permute(Mo(:,:,:) .* brainmask, [3 2 1]);
im = flip(im,1);
imshow(im(:,:,axial), [0 prctile(Mo(:), 98)])
axis square
colormap(hax(2, 1), 'gray')
ylabel('M0', 'rotation',0,'VerticalAlignment','middle');

axes(hax(2, 2));
im = permute(Mo(:,:,:) .* brainmask, [1 2 3]);
imshow(im(:,:,coronal), [0 prctile(Mo(:), 98)])
axis square
colormap(hax(2, 2), 'gray')

axes(hax(2, 3));
im = permute(Mo(:,:,:) .* brainmask, [1 3 2]);
im = flip(im, 2);
imshow(im(:,:,saggital), [0 prctile(Mo(:), 98)])
axis square
colormap(hax(2, 3), 'gray')
pos = get(hax(2, 3), 'Position');
hc = colorbar;
ylabel(hc, '(unitless)')
set(hax(2, 3), 'Position', pos)

% %%%%%%%%%%%%%% T1 (axial, sagittal, coronal)
load('GRrecon_T1Mo.mat', 'T1')
cm = load('/server/home/dce/STARDCE/T1estimation/Plotting_Functions/T1cm.mat', 'T1colormap');
% cm = load('/Users/yannickbliesener/Research/DCE/STARDCE/T1estimation/Plotting_Functions/T1cm.mat', 'T1colormap');
cm = cm.T1colormap;

axes(hax(3, 1));
im = permute(T1(:,:,:) / 1000 .* brainmask, [3 2 1]);
im = flip(im,1);
imshow(im(:,:,axial), [0 4.2])
axis square
colormap(hax(3, 1), cm)
ylabel('T1', 'rotation',0,'VerticalAlignment','middle');

axes(hax(3, 2));
im = permute(T1(:,:,:) / 1000 .* brainmask, [1 2 3]);
imshow(im(:,:,coronal), [0 4.2])
axis square
colormap(hax(3, 2), cm)

axes(hax(3, 3));
im = permute(T1(:,:,:) / 1000 .* brainmask, [1 3 2]);
im = flip(im, 2);
imshow(im(:,:,saggital), [0 4.2])
axis square
colormap(hax(3, 3), cm)
pos = get(hax(3, 3), 'Position');
hc = colorbar;
ylabel(hc, '[s]')
set(hax(3, 3), 'Position', pos)

% %%%%%%%%%%%%%% Pre (axial, sagittal, coronal)
load('GRrecon_perm_pre_post.mat', 'imgPRE')

axes(hax(4, 1));
im = permute(imgPRE(:,:,:) .* brainmask, [3 2 1]);
im = flip(im,1);
imshow(im(:,:,axial), [0 prctile(imgPRE(:), 98)])
axis square
colormap(hax(4, 1), 'gray')
ylabel('PreC', 'rotation',0,'VerticalAlignment','middle');

axes(hax(4, 2));
im = permute(imgPRE(:,:,:) .* brainmask, [1 2 3]);
imshow(im(:,:,coronal), [0 prctile(imgPRE(:), 98)])
axis square
colormap(hax(4, 2), 'gray')

axes(hax(4, 3));
im = permute(imgPRE(:,:,:) .* brainmask, [1 3 2]);
im = flip(im, 2);
imshow(im(:,:,saggital), [0 prctile(imgPRE(:), 98)])
axis square
colormap(hax(4, 3), 'gray')
pos = get(hax(4, 3), 'Position');
hc = colorbar;
ylabel(hc, '(unitless)')
set(hax(4, 3), 'Position', pos)

% %%%%%%%%%%%%%% Post (axial, sagittal, coronal)
load('GRrecon_perm_pre_post.mat', 'imgPOST')

axes(hax(5, 1));
im = permute(imgPOST(:,:,:) .* brainmask, [3 2 1]);
im = flip(im,1);
imshow(im(:,:,axial), [0 prctile(imgPOST(:), 98)])
axis square
colormap(hax(5, 1), 'gray')
ylabel('PostC', 'rotation',0,'VerticalAlignment','middle');

axes(hax(5, 2));
im = permute(imgPOST(:,:,:) .* brainmask, [1 2 3]);
imshow(im(:,:,coronal), [0 prctile(imgPOST(:), 98)])
axis square
colormap(hax(5, 2), 'gray')

axes(hax(5, 3));
im = permute(imgPOST(:,:,:) .* brainmask, [1 3 2]);
im = flip(im, 2);
imshow(im(:,:,saggital), [0 prctile(imgPOST(:), 98)])
axis square
colormap(hax(5, 3), 'gray')
pos = get(hax(5, 3), 'Position');
hc = colorbar;
ylabel(hc, '(unitless)')
set(hax(5, 3), 'Position', pos)


set(gcf, 'Position', [460 200 760 1140]);
set(gcf, 'Color', 'w')
set(findall(gcf, '-property','FontSize'),'FontSize',24)
export_fig(h2, 'spatial_maps.pdf', '-p0.02')
