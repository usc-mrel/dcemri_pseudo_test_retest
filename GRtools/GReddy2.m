function EC = GReddy2(kE,channels)

%   Number of points and receivers
np = size(kE,1);
nv = size(kE,2);
nr = size(kE,5);
if nargin < 2 || isempty(channels)
    channels = 1:nr;
end
r2d = 180/pi;

%   Strip unwanted channels
iE = kE(:,:,:,:,channels);
nr = length(channels);

%   Convert to image space in the readout direction
iE = iFastFT(fftshiftF(iE,1),1,1);

%   Split, extract, reshape, average
zE = iE(:,1:nv/2,:);
zE = reshape(zE,[np 3 2 8 24 nr]);
zE = squeeze(mean(zE(:,2,:,:,:,:),5));
yE = iE(:,(nv/2+1):nv,:);
yE = reshape(yE,[np 3 2 8 24 nr]);
yE = squeeze(mean(yE(:,2,:,:,:,:),5));
clear iE;

%   Separate into donor and baselines
zEb = squeeze(zE(:,2,:,:));
zE = squeeze(zE(:,1,:,:));
yEb = squeeze(yE(:,2,:,:));
yE = squeeze(yE(:,1,:,:));


%   Convert to y and z spatial location
yE = iFastFT(fftshiftF(yE,2),2,1);
yEb = iFastFT(fftshiftF(yEb,2),2,1);
zE = iFastFT(fftshiftF(zE,2),2,1);
zEb = iFastFT(fftshiftF(zEb,2),2,1);

%   Take complex difference with baseline, convert to angle
zEa = angle(zE.*conj(zEb));
yEa = angle(yE.*conj(yEb));

%   Setup surface polynomial fit
pts = np*8;
[x,y] = ndgrid(linspace(-0.5,0.5,np),linspace(-0.5,0.5,8));
X = [ones(pts,1) x(:) y(:)];
Wz = diag(abs(zE(:)+zEb(:)));
Wy = diag(abs(yE(:)+yEb(:)));

%   Fit z donor
E = pinv(X'*Wz*X)*X'*Wz;
pz = E*zEa(:);
zB0 = pz(1);
zGx = pz(2);
zGz = pz(3);

%   Fit y donor
E = pinv(X'*Wy*X)*X'*Wy;
py = E*yEa(:);
yB0 = py(1);
yGx = py(2);
yGy = py(3);

figure(1);
mesh(x,y,r2d*zEa);set(gca,'ZLim',[-30 30])
hold on;
surf(x,y,r2d*reshape(X*pz,[np 8]));set(gca,'ZLim',[-30 30])
figure(2);
mesh(x,y,r2d*yEa);set(gca,'ZLim',[-30 30])
hold on;
surf(x,y,r2d*reshape(X*py,[np 8]));set(gca,'ZLim',[-30 30])


%   Store as one output variable
EC.y.B0 = yB0;
EC.y.Gx = yGx;
EC.y.Gy= yGy;
EC.z.B0 = zB0;
EC.z.Gx = zGx;
EC.z.Gz= zGz;


end
