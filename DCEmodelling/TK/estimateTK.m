function [TK, predConc, OUT] = estimateTK(conc, AIF, initial_param, opt)
%function [TK, predConc, OUT] = estimateTK(conc, AIF, initial_param, opt)
%   general routine for TK estimation
%
%   Inputs:
%       conc            concentration time curves, [nx, ny, nz, nt]
%       AIF             arterial input function    [nt, 1]
%       inital_param    initial parameters for iterative estimators
%                       depending on the model:
%                           .vp
%                           .kt
%                           .ve
%                           .kep
%                           .fp
%                           .shift
%                       if [] this will use LLSQ for initial guess estimation
%       opt             options container
%       .size           [nx, ny, nz, nt]
%       .time           time vector
%       .AIF
%           .shifts     array of shifts that should be applied, decimal
%                       numbers lead to Akima spline interpolation
%       .TK
%           .models     cell array of models that should be fit, 
%                       e.g., {'null', 'vpmodel', 'patlak', 'tofts', 'etk', 'tcxm'}
%                       Note: if multiple are specified the most complex
%                       model has to be last!
%           .selectionCriterion
%                       a model selction criterion: 'SSE', 'AIC', 'BIC',
%                       'Balvays'
%           .parameterBounds
%                       parameter bounds
%   Output:
%       TK
%       	.vp         fractional blood plasma volume
%       	.kt         ktrans
%       	.ve         frational EES volume
%       	.kep        kep
%        	.fp         blood flow
%       	.shift      shift paramter, ie, the index of the entry in opt.AIF.shifts
%       predConc        estimated concentration time curves, [nx, ny, nz, nt]
%       OUT             output information container
%           .modelMask  mask of estimated model: index corresponding to opt.TK.models
%           .sse        estimation sum-of-squares
%
%
%   Example setup:
%       This setup can be generated with 
%           opt = TKESTIMATION_optset;     
%
%       or manually by adjusting the following example:
%
%         opt.TK.models               = {'patlak', 'tofts', 'etk', 'tcxm'};    % 'null' | 'vpmodel' | 'patlak' | 'etk' | 'tcxm'
%         opt.TK.selectionCriterion   = 'AIC';                                 % selection criterion for model selection: 'AIC' | 'BIC' | 'Balvays' | 'SSE'
%         opt.TK.parameterBounds.vp   = [0 1];                                 % [min max]
%         opt.TK.parameterBounds.ve   = [0.02 1];
%         opt.TK.parameterBounds.kt   = [0 inf];                                % in 1/s
%         opt.TK.parameterBounds.kep  = [0 inf];
%         opt.TK.parameterBounds.fp   = [0.001 100];
%         opt.TK.tofts.initialFit     = 'LLSQ';                                 % determines the initialization for non-linear fitting of tofts:
%                                                                               % 'patlak' | 'LLSQ'
%                                                                               % 'LLSQ' is linear fitting as described in Murase MRM 2004 and Flouri MRM 2016
%         opt.TK.tofts.solver         = 'gpufit';                               % fitting routine for Tofts model
%                                                                               % 'gpufit' | 'LLSQ' | 'ncg'
%
%         opt.TK.etk.initialFit       = 'LLSQ';                                 % determines the initialization for non-linear fitting of etk/tcxm:
%                                                                               % 'patlak' | 'LLSQ'
%                                                                               % 'LLSQ' is linear fitting as described in Murase MRM 2004 and Flouri MRM 2016
%         opt.TK.etk.solver           = 'gpufit';                               % fitting routine for ETK model
%                                                                               % 'gpufit' | 'LLSQ' | 'ncg' |'gaussnewton'
%         opt.TK.etk.init.kep         = [];                                     % [] | array of values for kep in 1/s
%                                                                               % if not empty overwrites the value from opt.TK.etk.initialFit for kep and
%                                                                               % refits to retain best fit
%                                                                               % Example:
%                                                                               % opt.TK.etk.init.kep = linspace(0, 1/60, 10);
%                                                     
%          opt.TK.tcxm.initialFit      = 'LLSQ';                                % determines the initialization for non-linear fitting of etk/tcxm:
%                                                                               % 'LLSQ'
%                                                                               % 'LLSQ' is linear fitting as described in Murase MRM 2004 and Flouri MRM 2016
%          opt.TK.tcxm.solver          = 'gpufit';                              % fitting routine for 2CXM model
%          opt.AIF.shifts              = [0 1 -1];                              % time bin shifts that are applied prior to fitting, best fit is chosen
%                                                                               %   postive shifts AIF right
%                                                                               %   negative shifts AIF left
%
%
%   See also TKESTIMATION_optset
%
% Yannick Bliesener 2019
%

N = prod(opt.size(1:3));
Nt = opt.size(4);

conc = reshape(conc, N, Nt);

if length(AIF) ~= Nt
    % this can happen if a population AIF is loaded
    error('AIF has length %i, while concentration time curves have length %i!', length(AIF), Nt) 
end
AIF = reshape(AIF, Nt, 1);

defaultTK.ve = 0.5;
% under uniform U(0,1) (agnostic) prior, 0.5 is the MMSE estimator of
% ve ;)
defaultTK.fp = inf;

% computed shifted AIFs
Nmodels = length(opt.TK.models);
Nshifts = 1;
if ~isempty(opt.AIF.shifts)
    Nshifts = length(opt.AIF.shifts);
    shiftedAIFs = zeros(Nt, Nshifts);
    for s=1:Nshifts
        dtime = opt.time(2)-opt.time(1);
        shiftedAIFs(:,s) = shiftAIF(AIF, opt.time, opt.AIF.shifts(s)*dtime);
    end
else
    shiftedAIFs = AIF;
end

TK = [];

% fit the different models
for model = opt.TK.models
    
    switch model{1}
        
        case 'null'
            fitTK.null.vp  = zeros(N,1);
            fitTK.null.ve  = defaultTK.ve*ones(size(fitTK.null.vp));
            fitTK.null.kt  = zeros(N,1);
            fitTK.null.kep = zeros(N,1);
            fitTK.null.fp  = defaultTK.fp*ones(N,1);
            fitTK.null.shift = zeros(N,1);
            
            predConc = zeros([N Nt]);
            SSE = sum((predConc - conc).^2, 2);
            
            num_param = 0;
            
        case 'vpmodel'
            
            SSE = -ones(N,Nshifts);
            
            fitTK.vpmodel.vp = zeros(N,Nshifts);
            
            for s=1:Nshifts
                [ fitTK.vpmodel.vp(:,s) ] = conc2tk_vp(N, Nt, conc, shiftedAIFs(:,s), 0 );
                
                fitTK.vpmodel = projectTKparamToBounds(fitTK.vpmodel, opt.TK.parameterBounds);
                
                intAIF = zeros(size(AIF));
                predConc = model_patlak(N, Nt, fitTK.vpmodel.vp(:,s), zeros(N,1), shiftedAIFs(:,s), intAIF);
                SSE(:,s) = sum((predConc - conc).^2, 2);
            end
            
            [SSE, fitTK.vpmodel.shift] = min(SSE, [], 2);
            
            ind = sub2ind([N Nshifts], (1:N)', fitTK.vpmodel.shift(:));
            fitTK.vpmodel.vp = fitTK.vpmodel.vp( ind(:) );
            
            fitTK.vpmodel.kt = zeros(size(fitTK.vpmodel.vp));
            fitTK.vpmodel.kep = zeros(size(fitTK.vpmodel.vp));
            fitTK.vpmodel.ve  = defaultTK.ve*ones(size(fitTK.vpmodel.vp));
            fitTK.vpmodel.fp  = defaultTK.fp*ones(size(fitTK.vpmodel.vp));
            
            num_param = 1;
            
        case 'patlak'
            
            SSE = -ones(N,Nshifts);
            
            fitTK.patlak.vp = zeros(N,Nshifts);
            fitTK.patlak.kt = zeros(N,Nshifts);
            
            for s=1:Nshifts
                modelParam.AIF = shiftedAIFs(:,s);
                modelParam.time = opt.time;
                modelParam.time = modelParam.time - modelParam.time(1);
                modelParam.Integration = opt.AIF.integration;
                
                intAIF = integrateAIF(modelParam.time, shiftedAIFs(:,s), modelParam.Integration);
                
                [ fitTK.patlak.kt(:,s), fitTK.patlak.vp(:,s) ] = conc2tk_patlak(N, Nt, conc, shiftedAIFs(:,s), intAIF, 0 );
                
                fitTK.patlak = projectTKparamToBounds(fitTK.patlak, opt.TK.parameterBounds);
                
                predConc = model_patlak(N, Nt, fitTK.patlak.vp(:,s), fitTK.patlak.kt(:,s), shiftedAIFs(:,s), intAIF);
                SSE(:,s) = sum((predConc - conc).^2, 2);
            end
            
            [SSE, fitTK.patlak.shift] = min(SSE, [], 2);
            
            ind = sub2ind([N Nshifts], (1:N)',fitTK.patlak.shift(:));
            fitTK.patlak.vp = fitTK.patlak.vp( ind(:) );
            fitTK.patlak.kt = fitTK.patlak.kt( ind(:) );
            
            fitTK.patlak.kep = zeros(size(fitTK.patlak.kt));
            fitTK.patlak.fp  = defaultTK.fp*ones(size(fitTK.patlak.kt));
            fitTK.patlak.ve  = defaultTK.ve*ones(size(fitTK.patlak.kt));
            
            num_param = 2;
            
        case 'tofts'
            SSE = -ones(N,Nshifts);
            
            fitTK.tofts.vp  = zeros(N,Nshifts);
            fitTK.tofts.kt  = zeros(N,Nshifts);
            fitTK.tofts.ve  = defaultTK.ve*ones(N,Nshifts);
            fitTK.tofts.kep = zeros(N,Nshifts);
            fitTK.tofts.fp  = defaultTK.fp*ones(N,1);
            
            for s=1:Nshifts
                modelParam.AIF = shiftedAIFs(:,s);
                modelParam.time = opt.time;
                modelParam.time = modelParam.time - modelParam.time(1);
                modelParam.Integration = opt.AIF.integration;
                
                init.kt  = zeros(N,1);
                init.kep = zeros(N,1);
                init.ve  = defaultTK.ve*ones(N,1);
                
                % generate initial guess
                if ~ismember(opt.TK.tofts.solver, {'LLSQ'})
                    % black list non-iterative fitting routines that do not
                    % require initialization
                    
                    if isempty(initial_param) && ~ismember(opt.TK.tofts.solver, {'LLSQ'})
                        
                        
                        switch opt.TK.tofts.initialFit
                            
                            case 'patlak'
                                intAIF = integrateAIF(modelParam.time, shiftedAIFs(:,s), modelParam.Integration);
                                init.kt  = reshape(intAIF \ conc', N, 1);
                                init.kep = zeros(N,1);
                                init.ve  = defaultTK.ve * ones(N,1);
                            case 'LLSQ'
                                
                                if isfield(opt, 'brainmask')
                                    ind = find( opt.brainmask(:) );
                                else
                                    ind = 1:N;
                                end
                                
                                fitOpt = getFitOptions('LLSQ', []);
                                [ tmp_initkt, tmp_initkep  ] = conc2tk_tofts(length(ind), Nt, conc(ind,:), modelParam, fitOpt);
                                
                                init.kt(ind)  = tmp_initkt;
                                init.kep(ind) = tmp_initkep;
                                init.ve(ind)  = tmp_initkt ./ tmp_initkep;
                                init.ve( isnan(init.ve) | isinf(init.ve) ) = defaultTK.ve;
                                clear tmp_init*
                                
                            otherwise
                                error('Unknown initial fitting method when trying to fit Tofts model')
                        end
                    else
                        init.kt = initial_param.kt;
                        init.ve = initial_param.ve;
                        init.kep = initial_param.kep;
                    end
                end
                
                init = projectTKparamToBounds(init, opt.TK.parameterBounds);
                
                switch opt.TK.etk.solver
                    case 'gpufit'
                        if isfield(opt, 'brainmask')
                            ind = find( opt.brainmask(:) );
                        else
                            ind = 1:N;
                        end
                        
                        fitOpt.init.kt = init.kt;
                        fitOpt.init.ve = init.ve;
                        fitOpt.init.kep = init.kep;
                        fitOpt.constraints = opt.TK.parameterBounds;
                        
                        for fn = fieldnames(fitOpt.constraints)'
                            fitOpt.constraints.(fn{1})(isinf(fitOpt.constraints.(fn{1}))) = 100;
                            % this line is important because the CUDA code
                            % cannot deal with inf ... so just set it something
                            % large
                        end
                        
                        fitOpt.init.kt = fitOpt.init.kt( ind(:) );
                        fitOpt.init.ve = fitOpt.init.ve( ind(:) );
                        fitOpt.init.kep = fitOpt.init.kep( ind(:) );
                        

                        [ tmp_kt, tmp_kep, flag, res, iter ] = conc2tk_tofts_gpu(length(ind(:)), Nt, conc(ind(:), :), modelParam, fitOpt);
                        
                        fitTK.tofts.kt(ind(:),s)  = tmp_kt;
                        fitTK.tofts.kep(ind(:),s) = tmp_kep;
                        clear tmp_kt tmp_kep
                        
                        fitTK.tofts.kt( isnan(fitTK.tofts.kt(:)) )   = 0;
                        fitTK.tofts.kep( isnan(fitTK.tofts.kep(:)) ) = 0;
                        
                        fitTK.tofts.ve(:,s) = fitTK.tofts.kt(:,s) ./ fitTK.tofts.kep(:,s);
                        fitTK.tofts.ve( isnan(fitTK.tofts.ve(:)) ) = defaultTK.ve;
                        
                        fitTK.tofts = projectTKparamToBounds(fitTK.tofts, opt.TK.parameterBounds);
                        
                    case 'LLSQ'
                        
                        fitTK.tofts.kt(:,s)  = init.kt;
                        fitTK.tofts.kep(:,s) = init.kep;
                        fitTK.tofts.ve(:,s)  = init.ve;
                        
                        if isfield(opt, 'brainmask')
                            ind = find( opt.brainmask(:) );
                        else
                            ind = 1:N;
                        end
                        
                        fitOpt = getFitOptions('LLSQ', []);
                        [ tmp_initkt, tmp_initkep  ] = conc2tk_tofts(length(ind), Nt, conc(ind,:), modelParam, fitOpt);
                        
                        fitTK.tofts.kt(ind(:),s)  = tmp_initkt;
                        fitTK.tofts.kep(ind(:),s) = tmp_initkep;
                        fitTK.tofts.ve(ind(:),s)  = tmp_initkt ./ tmp_initkep;
                        fitTK.tofts.ve( isnan(fitTK.tofts.ve) | isinf(fitTK.tofts.ve) ) = defaultTK.ve;
                        clear tmp_init*
                        
                        fitTK.tofts = projectTKparamToBounds(fitTK.tofts, opt.TK.parameterBounds);
                        
                        
                    otherwise
                        error('Unkown solver for Tofts model!')
                end
                
                predConc = model_standard(N, Nt, zeros(size(fitTK.tofts.kt(:,s))), fitTK.tofts.kt(:,s), fitTK.tofts.kep(:,s), shiftedAIFs(:,s), modelParam.time, modelParam.Integration);
                SSE(:,s) = sum((predConc - conc).^2, 2);
            end
            
            [SSE, fitTK.tofts.shift] = min(SSE, [], 2);
            
            % for each shifted AIF pick the best fit
            ind = sub2ind([N Nshifts], (1:N)',fitTK.tofts.shift(:));
            fitTK.tofts.kt  = fitTK.tofts.kt( ind(:) );
            fitTK.tofts.kep = fitTK.tofts.kep( ind(:) );
            fitTK.tofts.ve  = fitTK.tofts.ve( ind(:) );
            
            num_param = 2;
            
        case 'etk'
            
            SSE = -ones(N,Nshifts);
            
            fitTK.etk.vp  = zeros(N,Nshifts);
            fitTK.etk.kt  = zeros(N,Nshifts);
            fitTK.etk.ve  = defaultTK.ve*ones(N,Nshifts);
            fitTK.etk.kep = zeros(N,Nshifts);
            fitTK.etk.fp  = defaultTK.fp*ones(N,1);
            
            for s=1:Nshifts
                modelParam.AIF = shiftedAIFs(:,s);
                modelParam.time = opt.time;
                modelParam.time = modelParam.time - modelParam.time(1);
                modelParam.Integration = opt.AIF.integration;
                
                init.kt  = zeros(N,1);
                init.vp  = zeros(N,1);
                init.kep = zeros(N,1);
                init.ve  = defaultTK.ve*ones(N,1);
                
                % generate initial guess by fitting to patlak model
                if ~ismember(opt.TK.etk.solver, {'LLSQ'})
                    % black list non-iterative fitting routines that do not
                    % require initialization
                    
                    if isempty(initial_param) && ~ismember(opt.TK.etk.solver, {'LLSQ'})
                        
                        
                        switch opt.TK.etk.initialFit
                            
                            case 'patlak'
                                intAIF = integrateAIF(modelParam.time, shiftedAIFs(:,s), modelParam.Integration);
                                [ init.kt, init.vp ] = conc2tk_patlak(N, Nt, conc, shiftedAIFs(:,s), intAIF, 0 );
                                init.ve = defaultTK.ve*ones(N,1);
                                init.kep = init.kt ./ init.ve;
                            case 'LLSQ'
                                
                                if isfield(opt, 'brainmask')
                                    ind = find( opt.brainmask(:) );
                                else
                                    ind = 1:N;
                                end
                                
                                fitOpt = getFitOptions('LLSQ', []);
                                [ tmp_initkt, tmp_initvp, tmp_initkep  ] = conc2tk_standard(length(ind), Nt, conc(ind,:), modelParam, fitOpt);
                                
                                init.kt(ind)  = tmp_initkt;
                                init.vp(ind)  = tmp_initvp;
                                init.kep(ind) = tmp_initkep;
                                init.ve(ind)  = tmp_initkt ./ tmp_initkep;
                                init.ve( isnan(init.ve) | isinf(init.ve) ) = defaultTK.ve;
                                clear tmp_init*
                                
                            otherwise
                                error('Unknown initial fitting method when trying to fit ETK')
                        end
                    else
                        init.kt = initial_param.kt;
                        init.ve = initial_param.ve;
                        init.vp = initial_param.vp;
                        init.kep = initial_param.kep;
                    end
                end
                
                init = projectTKparamToBounds(init, opt.TK.parameterBounds);
                
                switch opt.TK.etk.solver
                    case 'gpufit'
                        if isfield(opt, 'brainmask')
                            ind = find( opt.brainmask(:) );
                        else
                            ind = 1:N;
                        end
                        
                        fitOpt.init.kt = init.kt;
                        fitOpt.init.ve = init.ve;
                        fitOpt.init.vp = init.vp;
                        fitOpt.init.kep = init.kep;
                        fitOpt.constraints = opt.TK.parameterBounds;
                        
                        for fn = fieldnames(fitOpt.constraints)'
                            fitOpt.constraints.(fn{1})(isinf(fitOpt.constraints.(fn{1}))) = 100;
                            % this line is important because the CUDA code
                            % cannot deal with inf ... so just set it something
                            % large
                        end
                        
                        fitOpt.init.kt = fitOpt.init.kt( ind(:) );
                        fitOpt.init.ve = fitOpt.init.ve( ind(:) );
                        fitOpt.init.vp = fitOpt.init.vp( ind(:) );
                        fitOpt.init.kep = fitOpt.init.kep( ind(:) );
                        
                        if 0
                            
                            [ tmp_kt, tmp_vp, tmp_ve, flag, res, iter ] = conc2tk_etofts_gpu(length(ind(:)), Nt, conc(ind(:), :), modelParam, fitOpt);
                            
                            fitTK.etk.kt(ind(:),s) = tmp_kt;
                            fitTK.etk.vp(ind(:),s) = tmp_vp;
                            fitTK.etk.ve(ind(:),s) = tmp_ve;
                            clear tmp_kt tmp_vp tmp_ve
                            
                            fitTK.etk.kt( isnan(fitTK.etk.kt(:)) ) = 0;
                            fitTK.etk.vp( isnan(fitTK.etk.vp(:)) ) = 0;
                            fitTK.etk.ve( isnan(fitTK.etk.ve(:)) ) = defaultTK.ve;
                            
                            fitTK.etk = projectTKparamToBounds(fitTK.etk, opt.TK.parameterBounds);
                            
                            fitTK.etk.kep(:,s) = fitTK.etk.kt(:,s) ./ fitTK.etk.ve(:,s);
                        else
                            
                            if isfield(opt.TK.etk, 'init') && isfield(opt.TK.etk.init, 'kep')  && ~isempty(opt.TK.etk.init.kep)
                                
                                prev_error = inf( size(ind) );
                                
                                tmp_kt  = fitOpt.init.kt;
                                tmp_vp  = fitOpt.init.vp;
                                tmp_kep = fitOpt.init.kep;
                                
                                for tmp_init_kep = opt.TK.etk.init.kep
                                    fitOpt.init.kep = tmp_init_kep * ones(size(fitOpt.init.vp));
                                    
                                    [ tmp_kt_, tmp_vp_, tmp_kep_, flag, res, iter ] = conc2tk_standard_gpu(length(ind(:)), Nt, conc(ind(:), :), modelParam, fitOpt);
                                    
                                    new_conc = model_standard(length(ind), Nt, tmp_vp_, tmp_kt_, tmp_kep_, shiftedAIFs(:,s), modelParam.time, modelParam.Integration);
                                    new_error = sum( abs(new_conc - conc(ind(:), :)).^2, 2);
                                    
                                    tmp_ind = find( new_error < prev_error );
                                    tmp_kt( tmp_ind )  =  tmp_kt_( tmp_ind );
                                    tmp_vp( tmp_ind )  =  tmp_vp_( tmp_ind );
                                    tmp_kep( tmp_ind ) =  tmp_kep_( tmp_ind );
                                end
                                clear tmp_kt_ tmp_vp_ tmp_kep_
                                
                            else
                                [ tmp_kt, tmp_vp, tmp_kep, flag, res, iter ] = conc2tk_standard_gpu(length(ind(:)), Nt, conc(ind(:), :), modelParam, fitOpt);
                            end
                            
                            fitTK.etk.kt(ind(:),s)  = tmp_kt;
                            fitTK.etk.vp(ind(:),s)  = tmp_vp;
                            fitTK.etk.kep(ind(:),s) = tmp_kep;
                            clear tmp_kt tmp_vp tmp_kep
                            
                            fitTK.etk.kt( isnan(fitTK.etk.kt(:)) )   = 0;
                            fitTK.etk.vp( isnan(fitTK.etk.vp(:)) )   = 0;
                            fitTK.etk.kep( isnan(fitTK.etk.kep(:)) ) = 0;
                            
                            fitTK.etk.ve(:,s) = fitTK.etk.kt(:,s) ./ fitTK.etk.kep(:,s);
                            fitTK.etk.ve( isnan(fitTK.etk.ve(:)) ) = defaultTK.ve;
                            
                            fitTK.etk = projectTKparamToBounds(fitTK.etk, opt.TK.parameterBounds);
                        end
                        
                    case 'LLSQ'
                        
                        fitTK.etk.kt(:,s)  = init.kt;
                        fitTK.etk.vp(:,s)  = init.vp;
                        fitTK.etk.kep(:,s) = init.kep;
                        fitTK.etk.ve(:,s)  = init.ve;
                        
                        if isfield(opt, 'brainmask')
                            ind = find( opt.brainmask(:) );
                        else
                            ind = 1:N;
                        end
                        
                        fitOpt = getFitOptions('LLSQ', []);
                        [ tmp_initkt, tmp_initvp, tmp_initkep  ] = conc2tk_standard(length(ind), Nt, conc(ind,:), modelParam, fitOpt);
                        
                        fitTK.etk.kt(ind(:),s)  = tmp_initkt;
                        fitTK.etk.vp(ind(:),s)  = tmp_initvp;
                        fitTK.etk.kep(ind(:),s) = tmp_initkep;
                        fitTK.etk.ve(ind(:),s)  = tmp_initkt ./ tmp_initkep;
                        fitTK.etk.ve( isnan(fitTK.etk.ve) | isinf(fitTK.etk.ve) ) = defaultTK.ve;
                        clear tmp_init*
                        
                        fitTK.etk = projectTKparamToBounds(fitTK.etk, opt.TK.parameterBounds);
                        
                    case {'ncg'}
                        kt = init.kt;
                        vp = init.vp;
                        kep = zeros([N 1]);
                        
                        if isfield(opt, 'brainmask')
                            ind = find( opt.brainmask(:) );
                        else
                            ind = 1:N;
                        end
                        
                        parfor n = ind
                            temp_fitOpt = getFitOptions('ncg', [vp(n); kt(n); 0])
                            [ kt(n), vp(n), kep(n)  ] = conc2tk_standard(1, Nt, conc(n,:), modelParam, temp_fitOpt);
                        end
                        fitTK.etk.kt(ind(:),s)  = kt;
                        fitTK.etk.vp(ind(:),s)  = vp;
                        fitTK.etk.kep(ind(:),s) = kep;
                        fitTK.etk.ve(ind(:),s)  = kt ./ kep;
                        fitTK.etk.ve( isnan(fitTK.etk.ve) | isinf(fitTK.etk.ve) ) = defaultTK.ve;
                        
                        fitTK.etk = projectTKparamToBounds(fitTK.etk, opt.TK.parameterBounds);
                        
                    case {'gaussnewton'}
                        fitTK.etk.kt(:,s)  = init.kt;
                        fitTK.etk.vp(:,s)  = init.vp;
                        fitTK.etk.kep(:,s) = init.kep;
                        fitTK.etk.ve(:,s)  = init.ve;
                        
                        if isfield(opt, 'brainmask')
                            ind = find( opt.brainmask(:) );
                        else
                            ind = 1:N;
                        end
                        
                        fitOpt = getFitOptions('gaussnewton', [init.vp(ind(:)); init.kt(ind(:)); init.kep(ind(:))]);
                        [ kt, vp, kep  ] = conc2tk_standard(length(ind), Nt, conc(ind,:), modelParam, fitOpt);
                        fitTK.etk.kt(ind(:),s)  = kt;
                        fitTK.etk.vp(ind(:),s)  = vp;
                        fitTK.etk.kep(ind(:),s) = kep;
                        fitTK.etk.ve(ind(:),s)  = kt ./ kep;
                        fitTK.etk.ve( isnan(fitTK.etk.ve) | isinf(fitTK.etk.ve) ) = defaultTK.ve;
                        
                        fitTK.etk = projectTKparamToBounds(fitTK.etk, opt.TK.parameterBounds);
                        
                    otherwise
                        error('Unkown solver for ETK model!')
                end
                
                predConc = model_standard(N, Nt, fitTK.etk.vp(:,s), fitTK.etk.kt(:,s), fitTK.etk.kep(:,s), shiftedAIFs(:,s), modelParam.time, modelParam.Integration);
                SSE(:,s) = sum((predConc - conc).^2, 2);
            end
            
            [SSE, fitTK.etk.shift] = min(SSE, [], 2);
            
            % for each shifted AIF pick the best fit
            ind = sub2ind([N Nshifts], (1:N)',fitTK.etk.shift(:));
            fitTK.etk.vp  = fitTK.etk.vp( ind(:) );
            fitTK.etk.kt  = fitTK.etk.kt( ind(:) );
            fitTK.etk.kep = fitTK.etk.kep( ind(:) );
            fitTK.etk.ve  = fitTK.etk.ve( ind(:) );
            
            num_param = 3;
            
            % this should not be necessary actually, yet as it turns out
            % the fits generated with Gpufit are rather poor sometimes and
            % can even lead to worsening of the fit
            % this is a hack to mend this
            if ~isempty(initial_param)
                prev_conc = tk2conc(initial_param, AIF, opt);
                prev_conc = reshape(prev_conc, N, Nt);
                prev_error = sum( abs(prev_conc - conc).^2, 2);
                
                new_conc = tk2conc(fitTK.etk, AIF, opt);
                new_conc = reshape(new_conc, N, Nt);
                new_error = sum( abs(new_conc - conc).^2, 2);
                
                ind = (new_error > prev_error);
                for fn = fieldnames(fitTK.etk)'
                    fitTK.etk.(fn{1})(ind(:)) = initial_param.(fn{1})(ind(:));
                end
                clear prev_conc prev_error new_conc new_error ind
            end
            
            
        case 'tcxm'
            % TODO: this is more like a draft
            warning('2CXM has never been tested before!')
            
            SSE = -ones(N,Nshifts);
            
            fitTK.tcxm.vp  = zeros(N,Nshifts);
            fitTK.tcxm.kt  = zeros(N,Nshifts);
            fitTK.tcxm.ve  = defaultTK.ve*ones(N,Nshifts);
            fitTK.tcxm.kep = zeros(N,Nshifts);
            fitTK.tcxm.fp  = defaultTK.fp*ones(N,Nshifts);
            
            for s=1:Nshifts
                modelParam.AIF = shiftedAIFs(:,s);
                modelParam.time = opt.time;
                modelParam.time = modelParam.time - modelParam.time(1);
                modelParam.Integration = opt.AIF.integration;
                
                init.kt  = zeros(N,1);
                init.vp  = zeros(N,1);
                init.fp  = defaultTK.fp*ones(N,1);
                init.ve  = defaultTK.ve*ones(N,1);
                
                if isfield(opt, 'brainmask')
                    ind = find( opt.brainmask(:) );
                else
                    ind = 1:N;
                end
                
                % generate initial guess by fitting
                if ~ismember(opt.TK.tcxm.solver, {'LLSQ'})
                    % black list non-iterative fitting routines that do not
                    % require initialization
                    
                    if isempty(initial_param)
                        % generate initial guess
                        
                        fitOpt = getFitOptions('LLSQ', []);
                        [ tmp_initkt, tmp_initvp, tmp_initkep, tmp_initfp] = conc2tk_2cxm(length(ind), Nt, conc(ind,:), modelParam, fitOpt);
                        init.kt(ind)  = tmp_initkt;
                        init.vp(ind)  = tmp_initvp;
                        init.kep(ind) = tmp_initkep;
                        init.fp(ind)  = tmp_initfp;
                        init.ve(ind)  = tmp_initkt ./ tmp_initkep;
                        init.ve( isnan(init.ve) | isinf(init.ve) ) = defaultTK.ve;
                        clear tmp_init*
                    else
                        init.kt  = initial_param.kt;
                        init.ve  = initial_param.ve;
                        init.vp  = initial_param.vp;
                        init.fp  = initial_param.fp;
                        init.kep = initial_param.kep;
                    end
                end
                
                init = projectTKparamToBounds(init, opt.TK.parameterBounds);
                
                switch opt.TK.tcxm.solver
                    case 'gpufit'
 
                        fitOpt.init.kt = init.kt;
                        fitOpt.init.ve = init.ve;
                        fitOpt.init.vp = init.vp;
                        fitOpt.init.fp = init.fp;
                        fitOpt.constraints = opt.TK.parameterBounds;
                        
                        for fn = fieldnames(fitOpt.constraints)'
                            fitOpt.constraints.(fn{1})(isinf(fitOpt.constraints.(fn{1}))) = 100;
                            % this line is important because the CUDA code
                            % cannot deal with inf ... so just set it something
                            % large
                        end
                        
                        fitOpt.init.kt = fitOpt.init.kt( ind(:) );
                        fitOpt.init.ve = fitOpt.init.ve( ind(:) );
                        fitOpt.init.vp = fitOpt.init.vp( ind(:) );
                        fitOpt.init.fp = fitOpt.init.fp( ind(:) );
                        
                        
                        [ tmp_kt, tmp_vp, tmp_ve, tmp_fp, flag, res, iter ] = conc2tk_2cxm_gpu(length(ind(:)), Nt, conc(ind(:), :), modelParam, fitOpt);
                        
                        fitTK.tcxm.kt(ind(:),s)  = tmp_kt;
                        fitTK.tcxm.vp(ind(:),s)  = tmp_vp;
                        fitTK.tcxm.ve(ind(:),s)  = tmp_ve;
                        fitTK.tcxm.fp(ind(:),s)  = tmp_fp;
                        clear tmp_kt tmp_vp tmp_kep tmp_fp
                        
                        fitTK.tcxm.kt( isnan(fitTK.tcxm.kt(:)) )   = 0;
                        fitTK.tcxm.vp( isnan(fitTK.tcxm.vp(:)) )   = 0;
                        fitTK.tcxm.ve( isnan(fitTK.tcxm.ve(:)) )   = defaultTK.ve;
                        fitTK.tcxm.fp( isnan(fitTK.tcxm.fp(:)) )   = defaultTK.fp;
                        
                        fitTK.tcxm.kep(:,s) = fitTK.tcxm.kt(:,s) ./ fitTK.tcxm.ve(:,s);
                        fitTK.tcxm.kep( isnan(fitTK.tcxm.kep(:)) ) = 0;
                        
                        fitTK.tcxm = projectTKparamToBounds(fitTK.tcxm, opt.TK.parameterBounds);
                        
                    case 'LLSQ'
                        fitTK.tcxm.kt(:,s)  = init.kt;
                        fitTK.tcxm.vp(:,s)  = init.vp;
                        fitTK.tcxm.fp(:,s)  = init.fp;
                        fitTK.tcxm.kep(:,s) = init.kep;
                        fitTK.tcxm.ve(:,s)  = init.ve;
                        
                        fitOpt = getFitOptions('LLSQ', []);
                        [ tmp_initkt, tmp_initvp, tmp_initkep, tmp_initfp  ] = conc2tk_2cxm(length(ind), Nt, conc(ind,:), modelParam, fitOpt);
                        
                        fitTK.tcxm.kt(ind(:),s)  = tmp_initkt;
                        fitTK.tcxm.vp(ind(:),s)  = tmp_initvp;
                        fitTK.tcxm.kep(ind(:),s) = tmp_initkep;
                        fitTK.tcxm.fp(ind(:),s)  = tmp_initfp;
                        fitTK.tcxm.ve(ind(:),s)  = tmp_initkt ./ tmp_initkep;
                        fitTK.tcxm.ve( isnan(fitTK.tcxm.ve) | isinf(fitTK.tcxm.ve) ) = defaultTK.ve;
                        clear tmp_init*
                        
                        fitTK.tcxm = projectTKparamToBounds(fitTK.tcxm, opt.TK.parameterBounds);
                        
                    otherwise
                        error('Unkown solver for 2CXM model!')
                end
                
                predConc = model_2cxm(N, Nt, fitTK.tcxm.vp(:,s), fitTK.tcxm.ve(:,s), fitTK.tcxm.kt(:,s), fitTK.tcxm.fp(:,s), shiftedAIFs(:,s), modelParam.time, modelParam.Integration);
                SSE(:,s) = sum((predConc - conc).^2, 2);
            end
            
            % pick the best shifted AIF among the results
            [SSE, fitTK.tcxm.shift] = min(SSE, [], 2);
            
            ind = sub2ind([N Nshifts], (1:N)',fitTK.tcxm.shift(:));
            fitTK.tcxm.vp  = fitTK.tcxm.vp( ind(:) );
            fitTK.tcxm.ve  = fitTK.tcxm.ve( ind(:) );
            fitTK.tcxm.kt  = fitTK.tcxm.kt( ind(:) );
            fitTK.tcxm.fp  = fitTK.tcxm.fp( ind(:) );
            fitTK.tcxm.kep = fitTK.tcxm.kep( ind(:) );
            
            num_param = 4;
            
        otherwise
            error('Unknown model to fit concentration time curves')
    end
    
    switch opt.TK.selectionCriterion
        
        % the 1 in all the ICs is added because the noise standard
        % deviation needs to be estimated too, hence the number of
        % free parameters is one higher than the TK model parameter
        case 'AIC'
            SC.(model{1}) = AIC_gaussian(num_param + 1, Nt, SSE, 1);
        case 'BIC'
            SC.(model{1}) = BIC_gaussian(num_param + 1, Nt, SSE, 1);
        case 'Balvays'
            % compute predicted concentration time curves
            predConc = tk2conc(fitTK.(model{1}), AIF, opt);
            
            [~, SC.(model{1})] = balvaysCriterion(conc,predConc);
            % FRI is best when its lowest and hence fits in the same
            % framework as AIC and BIC where lower is better!
        case 'SSE'
            SC.(model{1}) = SSE;
        otherwise
            error('Unkown model selection criterion')
    end
end

% select the best performing model for each voxel
if Nmodels == 1
    for fname = fieldnames(fitTK.(model{1}))'
        TK.(fname{1}) = fitTK.(model{1}).(fname{1});
    end
    mask = [];
else
    SC = struct2cell(SC);
    mask = minMask( SC{:});
    TKparam = fieldnames(fitTK.(model{1}))';
    for fname = TKparam
        TK.(fname{1}) = zeros([N 1]);
        for i = 1:Nmodels
            TK.(fname{1})(mask(:)==i) = fitTK.(opt.TK.models{i}).(fname{1})(mask(:)==i);
        end
    end
    mask = reshape(mask, opt.size(1:3));
end
OUT.modelMask = mask;

% compute output SSE
predConc = tk2conc(TK, AIF, opt);
predConc = reshape(predConc, N, Nt);
OUT.sse = sum((conc - predConc).^2, 2);
predConc = reshape(predConc, opt.size(1:4));

for fname = fieldnames(TK)'
    TK.(fname{1}) = reshape(TK.(fname{1}), opt.size(1:3));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%    Utilities  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function TK_ = projectTKparamToBounds(TK_, bounds)
        
        if ~isempty(bounds)
            for fname = fieldnames(TK_)'
                TK_.(fname{1}) = max(bounds.(fname{1})(1), TK_.(fname{1}));
                TK_.(fname{1}) = min(bounds.(fname{1})(2), TK_.(fname{1}));
            end
        end
        
    end


end

