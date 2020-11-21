function [k,petab,hdr,tstamp,pfle] = GRread2(path,mtrxL,mtrxH)

%   Use current directory as default
if nargin<1 || isempty(path)
    path = cd;
end
if nargin < 2
    mtrxL = [];
end
if nargin < 3
    mtrxH = [];
end
cd(path);

%   Obtain a directory listing
files = dir('P*.7');
name = files(1).name;

%   Open file
pfle = GERecon('Pfile.Load',name,'No-Anonymize');
hdr = GERecon('Pfile.Header');

%   Get some parameters
[np, nv, ns, nr, ph] = get_key_params2(pfle);
zip = hdr.RawHeader.zip_factor;
if zip > 1
    error('Zip factor > 1');
end

%   Read data
k = zeros([np nv ns nr ph],'single');
for phi = 1:ph
for nri = 1:nr
for nsi = 1:ns
    k(:,:,nsi,nri,phi) = GERecon('Pfile.KSpace',nsi,1,nri,phi);
end
end
end

%   Read data (option2)
%   the GERecon function seems to have an issue
%   The header size seems to be erroneous
%[~,k] = readGEdata(name);
k = double(k);

%   Reshape data
if numel(k) == (np*nv*ns*nr*ph)
    rvec1 = [np, nv, ns, nr, ph];
    pvec = [1, 3, 2, 5, 4];
    rvec2 = [np, nv*ns*ph, 1, 1, nr];
    k = reshape(k,rvec1);
elseif numel(k) == (np*(nv+1)*ns*nr*ph)
    rvec1 = [np, nv+1, ns, nr, ph];
    pvec = [1, 3, 2, 5, 4];
    rvec2 = [np, nv*ns*ph, 1, 1, nr];
    k = reshape(k,rvec1);
    k = k(:,2:nv+1,:,:,:);
else
    error('Unexpected data size');
end
k = flip(k,1);  %   This is totally empirical
k = permute(k,pvec);
k = reshape(k,rvec2);


%   Weight receiver coils based on noise
for i = 1:nr
    k(:,:,:,:,i) = k(:,:,:,:,i) ./ hdr.PrescanHeader.rec_std(i);
end


%   Read pe table
if exist('GRexttab.txt','file')
    petab = dlmread('GRexttab.txt');
else
    dirs = dir('*petable.txt');
    if length(dirs) == 1
        petab = dlmread(dirs.name);
    else
        error('Could not find PE table!')
    end
    
end
petab = petab(1:size(k,2),:);


%   Make timestamp
tr = hdr.ImageData.tr / 1e6;
tstamp = 0:tr:(size(k,2)*tr-tr);


%   Reduce resolution, if requested
if ~isempty(mtrxL)
    
    %   Check requested matrix size
    acqsize = [np nv ns];
    mtrxL = min(mtrxL,acqsize);
    mtrxL(mtrxL == 0) = acqsize(mtrxL == 0);
    
    %   Define new resolution
    xres = mtrxL(1);
    yres = mtrxL(2);
    zres = mtrxL(3);
    
    %   Crop readout
    k = k(round(np/2)-floor(xres/2)+1:round(np/2)+ceil(xres/2),:,1,1,:);
    
    %   Find phase encodes that match resolution requirements
    yl = round((nv-yres)/2)+1;
    yu = round((nv+yres)/2);
    zl = round((ns-zres)/2)+1;
    zu = round((ns+zres)/2);
    ind = find((petab(:,1) >= yl & petab(:,1) <= yu) & (petab(:,2) >= zl & petab(:,2)<=zu));
    
    %   Crop phase encodes
    k = k(:,ind,1,1,:);
    petab = petab(ind,:);
    petab(:,1) = petab(:,1) - round(nv/2) + round(yres/2);
    petab(:,2) = petab(:,2) - round(ns/2) + round(zres/2);
    tstamp = tstamp(ind);
    
    %   Update local variables
    np = xres;
    nv = yres;
    ns = zres;
    
    %   Update output variables
    pfle.xRes = np;
    pfle.yRes = nv;
    pfle.slices = ns;
    
end

%   Reconstruct on larger grid, if requested
if ~isempty(mtrxH)
    
    %   Check requested matrix size
    acqsize = [np nv ns];
    mtrxH = max(mtrxH,acqsize);
    mtrxH(mtrxH == 0) = acqsize(mtrxH == 0);
    
    xres = mtrxH(1);
    yres = mtrxH(2);
    zres = mtrxH(3);
    
    %   Pad readout
    if xres > np
        k2 = zeros([xres, size(petab,1), 1, 1, nr],class(k));
        indX = (round(xres/2) - round(np/2) + 1) : (round(xres/2) + round(np/2));
        k2(indX,:,:,:,:) = k;
        k = k2;
        clear k2;
    end
    
    %   Shift phase encode
    petab(:,1) = petab(:,1) - round(nv/2) + round(yres/2);
    petab(:,2) = petab(:,2) - round(ns/2) + round(zres/2);

    %   Update output variables
    pfle.xRes = xres;
    pfle.yRes = yres;
    pfle.slices = zres;
end

end

