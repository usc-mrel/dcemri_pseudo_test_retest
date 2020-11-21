function [tx,ty,tz,rx,ry,rz,img2] = estMOT(img,res,apply)

if nargin < 3
    apply = 0;
end
if nargin < 2
    res = [];
end
if nargin < 1
    error('Require at least one input');
end

%   Define some translation fit parameters
fitrotation = 0;
fittranslation = 1;
ups1 = 0.1;
ups2 = 0.01;
ups3 = 0.001;
mthd = 'spline';

%   Define some angular fit parameters
Nang  = 720  + 1;
Nangi = 3600 + 1;
Nrad  = 128;

%   Get size
[np, nv, ns, ph] = size(img);

%   Pad image in slice direction
img = cat(3,zeros([np nv 1 ph]),img,zeros([np nv 1 ph]));
ns = ns + 2;

%   Initialize motion array
tx = zeros(1,ph);ty = zeros(1,ph);tz = zeros(1,ph);
rx = zeros(1,ph);ry = zeros(1,ph);rz = zeros(1,ph);

if fittranslation
    %   Determine registration dimensions
    dim = [1 2 3];
    
    %   Convert to k-space
    k = fFastFT(img,dim,1);
    
    %   Get base image
    % kb = fFFT(mean(img,4),dim,1);
    kb = k(:,:,:,1);
    
    %   Precompute indices for translation fitting
    indx1 = (floor(np/2)+1-1):ups1:(floor(np/2)+1+1);
    indy1 = (floor(nv/2)+1-1):ups1:(floor(nv/2)+1+1);
    indz1 = (floor(ns/2)+1-1):ups1:(floor(ns/2)+1+1);
    indx2 = (floor(np/2)+1-1.1*ups1):ups2:(floor(np/2)+1+1.1*ups1);
    indy2 = (floor(nv/2)+1-1.1*ups1):ups2:(floor(nv/2)+1+1.1*ups1);
    indz2 = (floor(ns/2)+1-1.1*ups1):ups2:(floor(ns/2)+1+1.1*ups1);
    indx3 = (floor(np/2)+1-1.1*ups2):ups3:(floor(np/2)+1+1.1*ups2);
    indy3 = (floor(nv/2)+1-1.1*ups2):ups3:(floor(nv/2)+1+1.1*ups2);
    indz3 = (floor(ns/2)+1-1.1*ups2):ups3:(floor(ns/2)+1+1.1*ups2);
    sftx = floor(np/2)+1;
    sfty = floor(nv/2)+1;
    sftz = floor(ns/2)+1;
    [x,y,z] = ndgrid(1:np,1:nv,1:ns);
end


%   Baseline log-polar images for rotation
if fitrotation
    % ang = -pi:2*pi/(Nang-1):pi;
    % angi= -pi:2*pi/(Nangi-1):pi;
    % imgXb = squeeze(mean(img(round(np/2)-2:round(np/2)+2,:,:,1),1));
    % imgXb = imresize(imgXb,[512 512]);
    % imgXb = logpolar(imgXb,Nrad,Nang);
    % imgXb = fFFT(imgXb,[ 2],1);
    % imgYb = squeeze(mean(img(:,round(nv/2)-2:round(nv/2)+2,:,1),2));
    % imgYb = imresize(imgYb,[512 512]);
    % imgYb = logpolar(imgYb,Nrad,Nang);
    % imgYb = fFFT(imgYb,[ 2],1);
    % imgZb = squeeze(mean(img(:,:,round(ns/2)-2:round(ns/2)+2,1),3));
    % imgZb = imresize(imgZb,[512 512]);
    % imgZb = logpolar(imgZb,Nrad,Nang);
    % imgZb = fFFT(imgZb,[ 2],1);
end

%   Loop through time points
for int = 2:size(img,4)
    
    
    %%%% PERFORM TRANSLATION FIT %%%%
    if fittranslation
        %   Compute cross correlation
        C = kb.*conj(k(:,:,:,int));
        C = C./(abs(C)+eps);
        C = iFastFT(C,dim,0);
        C = fftshift(abs(C));
        
        %   Find peak on coarse index
        indC = find(C(:) == max(C(:)),1);
        [xC,yC,zC] = ind2sub([np nv ns],indC);
        xC = xC - sftx;
        yC = yC - sfty;
        zC = zC - sftz;
        
        
        %   Make fine grid and interpolate
        [X,Y,Z] = ndgrid(indx1+xC,indy1+yC,indz1+zC);
        Ci = interpn(x,y,z,C,X,Y,Z,mthd);
        
        %   Locate peak on fine grid
        ind = find(Ci(:) == max(Ci(:)),1);
        xC = X(ind) - sftx;
        yC = Y(ind) - sfty;
        zC = Z(ind) - sftz;
        
        
        %   Make superfine grid and interpolate
        [X,Y,Z] = ndgrid(indx2+xC,indy2+yC,indz2+zC);
        Ci = interpn(x,y,z,C,X,Y,Z,mthd);
        
        %   Locate peak on super fine grid
        ind = find(Ci(:) == max(Ci(:)),1);
        xC = X(ind) - sftx;
        yC = Y(ind) - sfty;
        zC = Z(ind) - sftz;
        
        
        %   Make super-super fine grid and interpolate
        [X,Y,Z] = ndgrid(indx3+xC,indy3+yC,indz3+zC);
        Ci = interpn(x,y,z,C,X,Y,Z,mthd);
        
        %   Locate peak on super-super fine grid
        ind = find(Ci(:) == max(Ci(:)),1);
        
        
        %   Save translation
        tx(int) = X(ind) - sftx;
        ty(int) = Y(ind) - sfty;
        tz(int) = Z(ind) - sftz;
    end
    %%%% DONE TRANSLATION %%%%
    
    
    
    %%%% PERFORM ROTATIONAL FITS %%%%
    if fitrotation
        %   X
        imgX = squeeze(mean(img(round(np/2)-2:round(np/2)+2,:,:,int),1));
        imgX = imresize(imgX,[512 512]);
        imgX = logpolar(imgX,Nrad,Nang);
        imgX = fFastFT(imgX,[ 2],1);
        C = imgXb.*conj(imgX);
        C = C./(abs(C)+eps);
        C = iFastFT(C,[ 2],1);
        C = fftshift(abs(C),2);
        C = sum(C,1);
        Ci = interp1(ang,C,angi,mthd);
        ind = find(Ci(:) == max(Ci(:)),1);
        rx(int) = angi(ind);
        %     ind = find(C(:) == max(C(:)),1);
        %     [~,t] = ind2sub([Nrad,Nang],ind);
        %     rx(int) = ang(t);
        
        %   Y
        imgY = squeeze(mean(img(:,round(nv/2)-2:round(nv/2)+2,:,int),2));
        imgY = imresize(imgY,[512 512]);
        imgY = logpolar(imgY,Nrad,Nang);
        imgY = fFastFT(imgY,[ 2],1);
        C = imgYb.*conj(imgY);
        C = C./(abs(C)+eps);
        C = iFFT(C,[ 2],1);
        C = fftshift(abs(C),2);
        C = sum(C,1);
        Ci = interp1(ang,C,angi,mthd);
        ind = find(Ci(:) == max(Ci(:)),1);
        ry(int) = angi(ind);
        %     ind = find(C(:) == max(C(:)),1);
        %     [~,t] = ind2sub([Nrad,Nang],ind);
        %     ry(int) = ang(t);
        
        %   Z
        imgZ = squeeze(mean(img(:,:,round(ns/2)-2:round(ns/2)+2,int),3));
        imgZ = imresize(imgZ,[512 512]);
        imgZ = logpolar(imgZ,Nrad,Nang);
        imgZ = fFastFT(imgZ,[ 2],1);
        C = imgZb.*conj(imgZ);
        C = C./(abs(C)+eps);
        C = iFastFT(C,[ 2],1);
        C = fftshift(abs(C),2);
        C = sum(C,1);
        Ci = interp1(ang,C,angi,mthd);
        ind = find(Ci(:) == max(Ci(:)),1);
        rz(int) = angi(ind);
        %     ind = find(C(:) == max(C(:)),1);
        %     [~,t] = ind2sub([Nrad,Nang],ind);
        %     rz(int) = ang(t);
    end
    
    
    %   Plot
    plot(1:int,tx(1:int),1:int,ty(1:int),1:int,tz(1:int),...
         1:int,180/pi*rx(1:int),'--',1:int,180/pi*ry(1:int),'--',1:int,180/pi*rz(1:int),'--');
    drawnow;
end

img2 = [];
if apply
    img2 = img;
    for int = 1:size(img,4)
        
        %   tx
        if tx(int) ~= 0
            ph = exp(-1i*pi*(-1:2/(np-1):1)*tx(int));
            ph = reshape(ph,[np 1 1]);
            ph = repmat(ph,[1 nv ns]);
            k(:,:,:,int) = k(:,:,:,int).*ph;
        end
        
        %   ty
        if ty(int) ~= 0
            ph = exp(-1i*pi*(-1:2/(nv-1):1)*ty(int));
            ph = reshape(ph,[1 nv 1]);
            ph = repmat(ph,[np 1 ns]);
            k(:,:,:,int) = k(:,:,:,int).*ph;
        end
        
        %   tz
        if tz(int) ~= 0
            ph = exp(-1i*pi*(-1:2/(ns-1):1)*tz(int));
            ph = reshape(ph,[1 1 ns]);
            ph = repmat(ph,[np nv 1]);
            k(:,:,:,int) = k(:,:,:,int).*ph;
        end
        
        %   Convert to image space
        img2(:,:,:,int) = iFastFT(k(:,:,:,int),dim,1);
        
        %   Get error measure
%         RMSEc(int) = sqrt(sum(sum(sum( abs(img(:,:,:,1) - img2(:,:,:,int)).^2))));
%         RMSE(int)  = sqrt(sum(sum(sum( abs(img(:,:,:,1) - img(:,:,:,int)).^2))));
        
    end
    
    img2 = img2(:,:,2:ns-1,:);
end

% ph = size(img,4);
% plot(1:ph,tx,1:ph,ty,1:ph,tz);
% figure;
% plot(1:ph,RMSE,1:ph,RMSEc);

%   Convert to spatial units
if nargin == 2 && ~isempty(res)
    rx = rx * res(1);
    ry = ry * res(2);
    rz = rz * res(3);
end

end


%%%% SUB-FUNCTIONS %%%%
function LPim = logpolar(im,nR,nT)

    [x,y] = size(im);
    rLim = min(y/2,x/2);
    rMin = rLim/4;
    Tran = linspace(0,pi*2,nT);
    Rran = exp(log(rMin):(log(rLim)-log(rMin))/(nR-1):log(rLim));
    [xRan,yRan] = meshgrid(Tran,Rran);
    xCrds = yRan.*cos(xRan)+x/2;
    yCrds = yRan.*sin(xRan)+y/2;
    
    [X,Y] = meshgrid((1:x),(1:y));
    LPim = interp2(X,Y,im,xCrds,yCrds,'bilinear');
    LPim(isnan(LPim)) = 0;

end



