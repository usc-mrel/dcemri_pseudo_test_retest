% Fermi filter
function H = FermiFilt(sz,N,T)

%   Dim 1
E = N(1)/sz(1);
x = -1:2/(sz(1)-1):1;
H1 = 1./(1+exp((abs(x)-E)/T))';
H = repmat(H1,[1,sz(2),sz(3)]);

%   Dim 2
E = N(2)/sz(2);
x = -1:2/(sz(2)-1):1;
H1 = reshape(1./(1+exp((abs(x)-E)/T)),[1 sz(2) 1]);
H = H.*repmat(H1,[sz(1),1,sz(3)]);

%   Dim3
E = N(3)/sz(3);
x = -1:2/(sz(3)-1):1;
H1 = reshape(1./(1+exp((abs(x)-E)/T)),[1 1 sz(3)]);
H = H.*repmat(H1,[sz(1),sz(2),1]);

end