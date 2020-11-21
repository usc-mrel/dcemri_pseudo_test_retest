%%
%   This script is the central reconstruction script of the STARDCE
%   package.
%
%
%
%
%   Yannick Bliesener, bliesene@Usc.edu
%   MREL at University of Southern California, 2020.
%

clear all; close all;

%% load dependencies
addpath( genpath('../orchestra-sdk-1.8-1.matlab') );
addpath( genpath('../dcemri_pseudo_test_retest') );
addpath( genpath('../Gpufit/Gpufit-build/Matlab') );

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
    
[T1, Mo, ~] = GRT1map(GRopt);

%% 2) reconstruct pre and post contrast (for brain mask)
GRopt = GRoptset(64, 'perm');

GRopt.mtrxH = matrixSize;

GRopt.DO_MOCCO = 0;
GRopt.DO_PRE_POST = 1;

[~, imgPRE, ~, imgPOST, ~, ~, ~] = GRrecon(GRopt, [], 'pre_post');

%% 3) reconstruct dynamic images
GRopt = GRoptset(64, 'perm');

GRopt.mtrxH = matrixSize;

GRopt.DO_MOCCO = 0;
GRopt.DO_SPSENSE = 1;

[~, ~, ~, ~, imgR, OUT] = GRrecon(GRopt, [], 'SPSENSE');

%% 4) estimate TK parameters

% sparse SENSE based TK estimation
[~, TK, TKESTOUT] = GRestimation(GRopt, [], 'SPSENSE');

% MOCCO based estimation
GRopt = GRoptset(64, 'perm');

GRopt.mtrxH = matrixSize;
    
[imgR, TK, OUT] = GRestimation(GRopt, path, 'MOCCO');

