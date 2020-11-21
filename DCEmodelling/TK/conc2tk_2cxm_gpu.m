function [ kt, vp, ve, fp, flag, res, iter ] = conc2tk_2cxm_gpu(N, Nt, conc, modelParam, opt)
%function [ kt, vp, ve, fp, flag, res, iter ] = conc2tk_2cxm_gpu(N, Nt, conc, modelParam, opt)
%   fits the PK parameters of the two-compartment exchange model
%
%   Inputs
%       N                   number of voxels
%       Nt                  number of time points
%       conc                concentration time curves [N, Nt]
%       modelParam
%           .time           time vector in s [Nt, 1]
%           .AIF            arterial input function [Nt, 1]
%       opt
%           .constraints    (optional) parameter constraints
%                           Don't put anything but a number for this! (ie no "inf")
%               .kt         [min, max] in 1/s
%               .ve         [min, max] 
%               .vp         [min, max]
%               .fp         [min, max] in ?ml/s/100 ml
%           .init           (optional) initial parameters
%               .kt         in 1/s [N, 1]
%               .ve         [N, 1]
%               .vp         [N, 1]
%               .fp         in ?ml/s/100 ml [N, 1]
%
%   Outputs
%       kt      in 1/s [N, 1]
%       vp      [N, 1]
%       ve      [N, 1]
%       fp      in ?ml/s/100 ml [N, 1]
%       flag
%           0:	The fit converged, tolerance is satisfied, the maximum number of iterations is not exceeded
%           1:	Maximum number of iterations exceeded
%           2:	During the Gauss-Jordan elimination the Hessian matrix is indicated as singular
%           3:	Non-positive curve values have been detected while using MLE (MLE requires only positive curve values)
%           4:	State not read from GPU Memory
%       res     Residual sum of squares
%       iter    Iterations
%
%   Yannick Bliesener 2019
%

num_parameters = 4;

% prepare data constainer
conc = reshape(conc, [N Nt]);
conc = conc';

% simple sanity checks to prevent mex hiccup
if ~isreal(conc)
    error('Concentration time curves have to be real for GPUfit!')
end
conc = real(conc);
conc = single(conc);

% prepare constraints
constraints = [];
if isfield(opt, 'parameterBounds')
    constraints = [opt.parameterBounds.kt(:) * 60; opt.parameterBounds.ve(:); opt.parameterBounds.vp(:); opt.parameterBounds.fp(:) * 60];
    constraints = repmat(constraints, [1 N]);
    constraints = single(constraints);
end

% prepare model id for extended Tofts model
model_id = ModelID.TWO_COMPARTMENT_EXCHANGE;

% initial parameters
if isfield(opt, 'init')
    initial_parameters = cat(1, reshape(opt.init.kt * 60, 1, N), reshape(opt.init.ve, 1, N), reshape(opt.init.vp, 1, N), reshape(opt.init.fp * 60, 1, N));
    initial_parameters = single(initial_parameters);
else
    initial_parameters = 0.1*ones([num_parameters N], 'single');
    % do not initialize to zero
end

% tolerance
tolerance = single(1e-6);

% max
max_n_iterations = int32(100);

% parameters to fit
parameters_to_fit = ones(num_parameters, 1, 'int32');

% estiamtor ID
estimator_id = EstimatorID.LSE;

% put time vector and AIF in user provided info
% convert time to minutes
user_info = single([modelParam.time(:) / 60; modelParam.AIF(:)]);


% fit
if isempty(constraints)
[parameters, flag, chi_squares, iter, time]...
    = gpufit(conc, [], model_id, initial_parameters, tolerance, max_n_iterations, parameters_to_fit, estimator_id, user_info);
else
[parameters, flag, chi_squares, iter, time]...
    = gpufit_constraints(conc, constraints, model_id, initial_parameters, tolerance, max_n_iterations, parameters_to_fit, estimator_id, user_info);
end
% extract the output
parameters = reshape(parameters, [num_parameters N]);
kt = reshape(parameters(1,:), N, 1) / 60;
ve = reshape(parameters(2,:), N, 1);
vp = reshape(parameters(3,:), N, 1);
fp = reshape(parameters(4,:), N, 1) / 60;

res = sum(chi_squares);

