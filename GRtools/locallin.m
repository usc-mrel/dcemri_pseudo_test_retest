function [yfit,slope,int,R2] = locallin(x,y,n)

%   Make sure n is odd
n = 2*floor(n/2)+1;

x = x(:);
y = y(:);

l = length(y);
fn2 = floor(n/2);
cn2 = ceil(n/2);

%   Initialize outputs
yfit = zeros(size(y),'like',y);
slope = zeros(size(y),'like',y);
int = zeros(size(y),'like',y);
ONE = ones([n,1],'like',x);

%   First points
xloc = x(1:n);
yloc = y(1:n);
A = [xloc ONE];
AI = pinv(A);
p = AI * yloc;
for i = 1:fn2
    xcur = x(i);
    yfit(i) = [xcur 1] * p;
    slope(i) = p(1);
    int(i) = p(2);
end

%   Middle points
for i = cn2:l-cn2
    ind = i-fn2:i+fn2;
    xloc = x(ind);
    yloc = y(ind);
    xcur = x(i);
    
    A = [xloc ONE];
    AI = pinv(A);
    p = AI * yloc;
    yfit(i) = [xcur 1] * p;
    slope(i) = p(1);
    int(i) = p(2);
end

%   Last points
xloc = x(l-n+1:l);
yloc = y(l-n+1:l);
A = [xloc ONE];
AI = pinv(A);
p = AI * yloc;
for i = l-fn2:l
    xcur = x(i);
    yfit(i) = [xcur 1] * p;
    slope(i) = p(1);
    int(i) = p(2);
end

%   Compute R2
SSres = sum(abs(y - yfit).^2);
SStot = sum(abs(y - mean(y)).^2);
R2 = 1-SSres/SStot;
