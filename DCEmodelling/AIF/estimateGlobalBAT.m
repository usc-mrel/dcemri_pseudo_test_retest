function tbolus = estimateGlobalBAT(dims, xtData, dtime, domain)
% function tbolus = estimateBAT(xtData)
%   estimates bolus arival time based on Marc's center of k-space method
%   
%   input:
%       xtData      image time series [Nx, Ny, Nz, Nt]
%       dtime       time difference between adjacent frames
%   output:
%       tbolus      BAT
%
%   Yannick 2018

Nx = dims(1);
Ny = dims(2);
Nz = dims(3);
Nt = dims(4);
N = Nx*Ny*Nz;

if nargin < 4
    domain = 'kspace';
end

switch domain
    case 'kspace'
        xtData = reshape(xtData, [Nx Ny Nz Nt]);
        ktData = fftnd( xtData, [1 2 3], 0);
        kCenter = abs(squeeze(ktData(1,1,1,:)) );

        [~,sl,in] = locallin(0:dtime:(Nt-1)*dtime, kCenter, 3);
        [~,tbolus] = findpeaks(sl, 'NPeaks', 1, 'MinPeakHeight', 0.5*max(sl));
    case 'imgspace'
        xtData = reshape(xtData, [N Nt]);
        tbolus = -ones(N,1);
        for i=1:N
            [~,sl,in] = locallin(0:dtime:(Nt-1)*dtime, xtData(i,:), 1);
            [~,tb] = findpeaks(sl, 'NPeaks', 1, 'MinPeakHeight', 0.5*max(sl));
            if ~isempty( tb )
                tbolus(i) = tb;
            end
        end
    otherwise
        error('Unkown domain')
end

end

