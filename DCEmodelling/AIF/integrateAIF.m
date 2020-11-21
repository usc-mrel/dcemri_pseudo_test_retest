function intAIF = integrateAIF(tModel, AIF, integration)
%function intAIF = integrateAIF(tModel, AIF, integration)
%   integrates any array of time intensity curves along the last
%   dimension
%
%
%
%   Yannick 2019

if nargin < 3
    integration = 'sum';
end

% check if time bins are uniform, the last one may contain rounding errors
uniform_time = (length(uniquetol(tModel(2:end-1) - tModel(1:end-2),0.001)) == 1);
deltat = max(tModel(2:end-1) - tModel(1:end-2));

sizeAIF = size(AIF);
Nt = length( tModel );
AIF = reshape(AIF, [], Nt);

intAIF = zeros(size(AIF));

switch ( integration )
    case 'trapz'
        for tt=1:Nt
            if tt == 1
                intAIF(:, tt) = 0;
            else
                intAIF(:, tt) = trapz( tModel( 1:tt ), AIF(:, 1:tt ), 2 );
            end
        end
    case 'sum'
        if uniform_time
            intAIF(:, 2:end) = deltat .* cumsum(AIF(:,1:end-1), 2);
        else
            for tt=1:Nt
                if tt == 1
                    intAIF(:, tt) = 0;
                else
                    intAIF(:, tt) = intAIF(:, tt-1) + (tModel( tt ) - tModel( tt-1 )) .* AIF(:, tt );
                end
            end
        end
    otherwise
        error('Unknown integration')
end

intAIF = reshape(intAIF, sizeAIF);

end

