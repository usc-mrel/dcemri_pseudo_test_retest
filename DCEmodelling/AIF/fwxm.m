function dt = fwxm(AIF, time, frac)
% dt = fwxm(AIF, time, prct)
%   compute the Full-Width-of-frac-Maximum
%
%   Yannick 2020
%


z = AIF > frac*max(AIF);
zrise = find(z,1,'first');
if isempty(zrise)
    zrise = 1;
end
zfall = find((~z)&((1:length(time))>zrise),1,'first');
if isempty(zfall)
    zfall = length(time);
end
dt = time(zfall) - time(zrise);

end

