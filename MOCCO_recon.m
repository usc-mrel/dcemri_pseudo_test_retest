function [img,outputInfo] = MOCCO_recon(k,predSignal,init,opt)


% Yannick Bliesener
% 10/21/2017

% default values
tol     = opt.tol;
maxIter = opt.MaxIter;
solver  = opt.solver;

Nx = opt.size(1);
Ny = opt.size(2);
Nz = opt.size(3);
Nt = opt.size(4);
Nc = opt.size(5);
Nk = numel(k);
opt.Nk = Nk;

if strcmp( solver, 'pcg' )
    %% Prep RHS
    b = backward_model([k(:); sqrt(opt.lambda)*predSignal(:)], opt );
    
    %% Solve
    [x, flag, res, iter] = pcg(@(x) mulA(x, opt), b, tol , maxIter, [], [], init(:) );
elseif strcmp( solver, 'lsqr' )
    %% Prep RHS
    b = [k(:); sqrt(lambda)*predSignal(:)];
    
    %% Solve
    [x, flag, res, iter] = lsqr(@(x, transA) mulA(x, opt, transA ) , b(:), tol , maxIter, [], [], init(:));
else
    error('recon_sense: Solver unknown!');
end

img = reshape( x, [Nx Ny Nz Nt]);

outputInfo.flag = flag;
outputInfo.res  = res;
outputInfo.iter = iter;


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
            y = forward_backward_model(x, opt);
        end
    end

    function [ kc ] = forward_model(x, opt)
        
        x = reshape(x, opt.size(1:4));
        
        % Convert to k-space
        kc = fSENSE(x,opt);
        
        switch opt.SENSE.kspaceDataType
            case 'full'
                kc = opt.fFT(kc,opt);
            case 'list'
                kc = opt.fFTU(kc,opt);
            otherwise
                error('Unkown k-space data format!')
        end
        

        %% Data consistency && MOCCO part
        kc = [kc(:); sqrt(opt.lambda)*x(:)];

    end

    function [ x ] = backward_model(y, opt )

        Nx = opt.size(1);
        Ny = opt.size(2);
        Nz = opt.size(3);
        Nt = opt.size(4);
        Nc = opt.size(5);

        % DC part
        kc = y(1:opt.Nk);
        
        switch opt.SENSE.kspaceDataType
            case 'full'
                kc = reshape(kc, opt.size);
                kc = opt.iFT(kc,opt);
            case 'list'
                kc = opt.iFTU(kc,opt);
            otherwise
                error('Unkown k-space data format!')
        end
        kc = kc * prod(opt.size(opt.FTdim));
        kc = opt.iPI(kc,opt);


        %% MOCCO part
        xm = y(opt.Nk+1:end);

        %% Vectorize
        x = kc(:) + sqrt(opt.lambda)*xm(:);
    end

    function [ kc ] = forward_backward_model(x, opt)
        
        x = reshape(x, opt.size(1:4));
        
        % Convert to k-space
        kc = fSENSE(x,opt);
        
        switch opt.SENSE.kspaceDataType
            case 'full'
                kc = opt.fFT(kc,opt);
                kc = reshape(kc, opt.size);
                kc = opt.iFT(kc,opt);
            case 'list'
                kc = opt.fFTU(kc,opt);
                kc = opt.iFTU(kc,opt);
            otherwise
                error('Unkown k-space data format!')
        end
        kc = kc * prod(opt.size(opt.FTdim));
        kc = opt.iPI(kc,opt);
        
        %% Data consistency && MOCCO part
        kc = kc(:) + opt.lambda*x(:);
        
    end


end

