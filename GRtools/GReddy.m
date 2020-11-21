function EC = GReddy(kE,channels)

%   Number of points and receivers
np = size(kE,1);
nr = size(kE,5);
if nargin < 2 || isempty(channels)
    channels = 1:nr;
end
r2d = 180/pi;

%   Convert to image space in the readout direction
iE = iFastFT(fftshiftF(kE,1),1,1);

%   Split and reshape
zE = iE(:,1:100,:);
zE = reshape(zE,[np 10 10 nr]);
yE = iE(:,101:200,:);
yE = reshape(yE,[np 10 10 nr]);
clear iE;

%   Remove first pass and the encoding step
yE = yE(:,2:10,2:10,:);
zE = zE(:,2:10,2:10,:);

%   Get baseline (last 4 points)
bl = cat(3,zE(:,6:9,:,:),yE(:,6:9,:,:));
zE = zE(:,1:5,:,:);
yE = yE(:,1:5,:,:);

%   Average eddy current trials
zE = squeeze(mean(zE,3));
yE = squeeze(mean(yE,3));

%   Average baselines, replicate to match size of z/y eddy array
bl = squeeze(mean(mean(bl,2),3));
bl = reshape(bl,[np 1 nr]);
bl = repmat(bl,[1 size(zE,2) 1]);


%   Use only requested channels
%   Useful for coil compression, which makes the first channel very good, others not
yE = yE(:,:,channels);
zE = zE(:,:,channels);
bl = bl(:,:,channels);
nr = length(channels);

%   Take complex difference with baseline, convert to angle
zEa = angle(zE.*conj(bl));
yEa = angle(yE.*conj(bl));

%   Define weighting vector
W = abs(sum(yE,2) + sum(zE,2) + sum(bl,2));
W = reshape(W,[np*nr 1]);
W = diag(W);

%   Fit a quadratic for each time point
x = linspace(0,1,np)';
X = [ones(np*nr,1) repmat(x,[nr 1]) repmat(x.^2,[nr 1])];
Xplot = [ones(np,1) x x.^2];
E = pinv(X'*W*X)*X'*W;
for i = 1:size(zEa,2)
    zEatmp = reshape(zEa(:,i,:),[np*nr 1]);
    pz = E*zEatmp;
    zB0(i) = pz(1);
    zGx(i) = pz(2);
    zGx2(i) = pz(3);
    
    yEatmp = reshape(yEa(:,i,:),[np*nr 1]);
    py = E*yEatmp;
    yB0(i) = py(1);
    yGx(i) = py(2);
    yGx2(i) = py(3);
    
    figure(1);
    plot(x,r2d*squeeze(zEa(:,i,:)),'b',...
         x,r2d*Xplot*pz,'b',...
         x,r2d*squeeze(yEa(:,i,:)),'r',...
         x,r2d*Xplot*py,'r');
    set(gca,'YLim',r2d*[-0.4 0.4])
    pause(0.1);
end

%   Fit eddy currents to exponential decay
tr = 0:(size(zEa,2)-1);
[yB0f.A, yB0f.t] = fit_eddy(tr,yB0);
[yGxf.A, yGxf.t] = fit_eddy(tr,yGx);
[yGx2f.A, yGx2f.t] = fit_eddy(tr,yGx2);
[zB0f.A, zB0f.t] = fit_eddy(tr,zB0);
[zGxf.A, zGxf.t] = fit_eddy(tr,zGx);
[zGx2f.A, zGx2f.t] = fit_eddy(tr,zGx2);

%   Store as one output variable
EC.y.B0 = yB0f;
EC.y.Gx = yGxf;
EC.y.Gx2= yGx2f;
EC.z.B0 = zB0f;
EC.z.Gx = zGxf;
EC.z.Gx2= zGx2f;

%   Plot
figure(2);
plot(tr,r2d*yB0,'bx',tr,r2d*yB0f.A*exp(-tr/yB0f.t),'b-',...
     tr,r2d*yGx,'rx',tr,r2d*yGxf.A*exp(-tr/yGxf.t),'r-',...
     tr,r2d*yGx2,'kx',tr,r2d*yGx2f.A*exp(-tr/yGx2f.t),'k-');
legend('yB0','yB0fit','yGx','yGxfit','yGx2','yGx2fit')
set(gca,'YLim',r2d*[-0.4 0.4]);
figure(3);
plot(tr,r2d*zB0,'bx',tr,r2d*zB0f.A*exp(-tr/zB0f.t),'b-',...
     tr,r2d*zGx,'rx',tr,r2d*zGxf.A*exp(-tr/zGxf.t),'r-',...
     tr,r2d*zGx2,'kx',tr,r2d*zGx2f.A*exp(-tr/zGx2f.t),'k-');
legend('zB0','zB0fit','zGx','zGxfit','zGx2','zGx2fit');
set(gca,'YLim',r2d*[-0.4 0.4]);

end



function [A,tau] = fit_eddy(tr,eddy)

%   Fit exponential decay
xo =  [  0   1.00];
xub = [ 10  10.00];
xlb = [-10   0.01];
options = optimset('TolFun',1e-6,'TolCon',1e-6,'MaxIter',10000,'MaxFunEvals',10000,...
                   'TolX',1e-6,'LargeScale','off','Display','off');

X = fmincon(@expobj,xo,[],[],[],[],xlb,xub,[],options,tr,double(eddy));

A = X(1);
tau = X(2);

function f = expobj(X,T,S)
    f = (sum((X(1).*exp(-T(:)./X(2)) - S(:)).^2));
end

end