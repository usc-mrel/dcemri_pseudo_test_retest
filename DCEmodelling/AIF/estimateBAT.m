function tbolus = estimateBAT(dims, conc, dtime, method)
% function tbolus = estimateBAT(xtData)
%   estimates bolus arival time based with an adaptation of Choeng et al.
%   
%   input:
%       conc        image time series [Nx, Ny, Nz, Nt]
%       dtime       time difference between adjacent frames
%       method      well, the method:
%                   - 'LL': piecewise linear approximation
%                   - 'LQ': piecewise linear-quadratic approximation
%   output:
%       tbolus      BAT
%
%   References:
%       [1] Cheong LH, Koh TS, Hou Z. An automatic approach for estimating bolus arrival time in dynamic contrast MRI using piecewise continuous 
%           regression models. Phys Med Biol. 2003;83(5):N83--N88. http://europepmc.org/abstract/MED/12696805.
%
%   Yannick 2018

Nx = dims(1);
Ny = dims(2);
Nz = dims(3);
Nt = dims(4);
N = Nx*Ny*Nz;

if nargin < 4
    method = 'LL';
end

globalBAT = estimateGlobalBAT(dims, conc, dtime, 'kspace');

conc = reshape(conc, N, Nt);
[~, initialPeakEstimate] = max(conc(:,1:globalBAT+5), [], 2); % the 5 is kind of arbitrary

tbolus = -ones([Nx Ny Nz]);

options = optimoptions('lsqlin','Display','off');

% for lsqlin
conc = double(conc);
dtime = double(dtime);

for i=1:N
    p = initialPeakEstimate(i);
    
    c = conc(i,:)';
    c = c(1:p);

    RMSE = zeros(p,1);
    for k = 1:p
        switch method
            case 'LL'
                % piecewise linear approximation
                X = [ones(p,1), zeros(p, 1)];
                X(k+1:end,2) = (1:p-k)*dtime;
                
                correction = 0; % fishy but works better
            case 'LQ'
                % piecewise linear-quadratic approximation
                X = [ones(p,1), zeros(p, 1), zeros(p, 1)];
                X(k+1:end,2) = (1:p-k)*dtime;
                X(k+1:end,3) = ((1:p-k)*dtime).^2;
                
                correction = 1; % fishy but works better
            otherwise
                error('estimateBAT: Unkown method')
        end
        
        [~,RMSE(k)] = lsqlin(X,c,[],[],[],[],zeros([size(X,2), 1]), [], zeros(size(X,2),1), options);
    end
    [~,tbolus(i)] = min(RMSE);
    tbolus(i) = tbolus(i) + correction;
end

end

