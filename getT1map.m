function [T1,Mo,b1r] = getT1map(sz,nalpha,getb1) %#ok<*AGROW>
%   Obtain B1 maps (dual angle) and T1 maps (multi flip angle)
%   
%   Author: RML
%   Date: 03/2014
%   
%   Usage: [T1,Mo,b1r] = getT1map(sz,nalpha,getb1)
%   
%   Input:
%   sz:     Desired image size. Match this to the resolution of the dynamic data
%   nalpha: Number of flip angles in T1 mapping
%   getb1:  Flag to also do B1+ mapping
%   
%   Output:
%   T1: T1 map (s)
%   Mo: Relative Mo maps (arb. units)
%   b1r: Relative B1+ maps (scale factor)


if nargin < 3
    getb1 = 1;
end
if nargin < 2
    nalpha = 3;
end
if nargin < 1
    error('Specify desired image size');
end

%   Get current directory
cdir = pwd;

%   Get directories
for i = 1:nalpha
    pname{i} = uigetdir(pwd,'Select T1 mapping directories');
    cd(pname{i});
    cd ..;
end
if getb1 == 1
    for i = 1:2
        pnameb1{i} = uigetdir(pwd,'Select B1 mapping images (alpha first, then 2*alpha)');
    end
elseif getb1 == 2
    pnameb1{1} = uigetdir(pwd,'Select B1 mapping image directory');
end

%   Get images
for i = 1:nalpha    
    cd(pname{i}); [par,a] = dicomr;
    cd ..;
    fa(i) = par.FlipAngle;
    tr(i) = par.RepetitionTime;
    img(:,:,:,i) = a;
end

%   Sort images by flip angle
[fa,ind] = sort(fa);
img = img(:,:,:,ind);
tr = tr(ind);
if any(tr ~= mean(tr))
    warning('TR is not constant: %g, %g, %g',tr(1),tr(2),tr(3));
end
tr = mean(tr);

%   Get relative B1 map
if getb1
    b1r = getB1map(pnameb1);
else
    b1r = ones([size(img,1) size(img,2) size(img,3)]);
end

%   Convert images to the desired size
for i = 1:nalpha
    img2(:,:,:,i) = imresize3d(img(:,:,:,i),sz);
end
b1r = imresize3d(b1r,sz);
clear img

%   Register volumes
% for i = 2:nalpha
%     [img2(:,:,:,i),~,motion(:,i)] = registervolume(img2(:,:,:,i),img2(:,:,:,1),...
%         [par.ReconstructionDiameter par.ReconstructionDiameter par.SliceThickness*size(img2,3)]);
% end
% plot(fa,motion);drawnow;

%   Make T1 map
[T1,Mo] = despot1(img2,fa,tr,b1r);

%   Filter to reduce edge effects
filt = gauss2(ones(7,7,5),0.75);
filt = filt./sum(filt(:));
wt = sum(img2,4);
T1t = convn(T1.*wt,filt,'same')./(convn(wt,filt,'same')+eps);
Mot = convn(Mo.*wt,filt,'same')./(convn(wt,filt,'same')+eps);
ind = find(T1>0);
T1t(ind) = T1(ind);
T1 = T1t;clear T1t wt;
Mot(ind) = Mo(ind);
Mo = Mot;clear Mot;

cd(cdir);
