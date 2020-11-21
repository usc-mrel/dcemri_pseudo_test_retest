function kcal = extract_kcal(k,n)

%   Get size
[np nv ns nt nr] = size(k);

%   Check n
if nargin < 2 || isempty(n)
    n = [np nv ns];
end

%   Extract the calibration data
if length(n) == 1
    n = [n nv ns];
elseif length(n) == 2
    n = [n(1) n(2) ns];
elseif length(n) > 3
    n = n(1:3);
end
indnp = round(np/2-n(1)/2):round(np/2+n(1)/2 - 1);
indpe = round(nv/2-n(2)/2):round(nv/2+n(2)/2 - 1);
indns = round(ns/2-n(3)/2):round(ns/2+n(3)/2 - 1);
kcal = mean(k(indnp,indpe,indns,1:nt,:),4);
