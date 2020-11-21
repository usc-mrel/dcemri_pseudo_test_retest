function AIF = estimateAIFfromTissueConcentration(conc, opt)
%function AIF = estimateAIFfromTissueConcentration(conc, opt)
%   fits AIF to concentration time curves given TK parameters
%
%   see estimateAIF.m
%
%   Yannick Bliesener, bliesene@usc.edu, 2019

% default values
tol     = opt.tol;
maxIter = opt.maxIter;
solver  = opt.solver;

initAIF = opt.init;

N = prod(opt.size(1:3));
Nt = opt.size(4);

if strcmp( solver, 'pcg' )
    %% Prep RHS
    b = backward_model(conc(:), opt);
    
    %% Solve
    [AIF, flag, res, iter] = pcg(@(x) mulA(x, opt), b, tol , maxIter, [], [], initAIF(:) );
elseif strcmp( solver, 'lsqr' )
    %% Prep RHS
    b = conc(:);
    
    %% Solve
    [AIF, flag, res, iter] = lsqr(@(x, transA) mulA(x, opt, transA ) , b, tol , maxIter, [], [], initAIF(:));
else
    error('recon_sense: Solver unknown!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Utilities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function y = mulA(x, opt, transA )
        if nargin > 2
            if strcmp( transA, 'notransp')
                y = forward_model(x, opt);
            elseif strcmp( transA, 'transp')
                y = backward_model(x, opt);
            else
                error('Unknown operator');
            end
        else
            z = forward_model(x, opt);
            y = backward_model(z, opt);
        end
    end

    function [ C ] = forward_model(AIF, opt)
        
        N = prod(opt.size(1:3));
        Nt = opt.size(4);

        vp = reshape(opt.TK.vp, N, 1);
        Ktrans = reshape(opt.TK.kt, N, 1);
        Kep = reshape(opt.TK.kep, N, 1);
        t = reshape(opt.time, [1 Nt]);
        
        AIF = reshape(AIF, [1 Nt]);

        C = vp*AIF;

        deltat = t(2) - t(1);
        F = exp( repmat(-Kep(:), [1 Nt]) .* repmat(t, [N 1]));
        
        AIF_fft = fft( [AIF, zeros(1,Nt-1)], [], 2);
        IRF_fft = fft( [F, zeros(N, Nt-1)], [], 2);
        Ct_ifft = ifft( repmat(AIF_fft, [N, 1]) .* IRF_fft, [], 2 );
        C = C + repmat( (deltat.*Ktrans(:)), [1 Nt]) .* Ct_ifft(:, 1:Nt);
        
        C = C(:);
        
    end

    function [ AIF ] = backward_model(C, opt )
        
        N = prod(opt.size(1:3));
        Nt = opt.size(4);
        
        vp = reshape(opt.TK.vp, N, 1);
        Ktrans = reshape(opt.TK.kt, N, 1);
        Kep = reshape(opt.TK.kep, N, 1);
        t = reshape(opt.time, [1 Nt]);
        deltat = t(2) - t(1);
        
        C = reshape(C, [N Nt]);

        AIF = repmat(vp, [1 Nt]) .* C;

        F = exp( repmat(-Kep(:), [1 Nt]) .* repmat(t, [N 1]));
        C_fft = fft( [C, zeros(N,Nt-1)], [], 2);
        IRF_fft = fft( [F, zeros(N, Nt-1)], [], 2);
        Ct_ifft = ifft( conj(C_fft) .* IRF_fft, [], 2 );
        AIF = AIF + repmat( (deltat.*Ktrans(:)), [1 Nt]) .* Ct_ifft(:, 1:Nt);
        AIF = sum(AIF, 1);
        
        AIF = AIF(:);
        
    end

end

