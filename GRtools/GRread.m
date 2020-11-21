function [k,petab,par,tstamp] = GRread(path,mtrxL,mtrxH,cls)

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
if nargin < 4
    cls = 'single';
end
cd(path);

%   Obtain a directory listing
files = dir('P*.7');
nfiles = length(files);

%   Loop through files
k = [];
for i = 1:nfiles
    
    name = files(i).name;
    
    %   Read data
    [par,ktmp] = readGEdata(name);
    if strcmp(cls,'single')
        ktmp = single(ktmp);
    else
        ktmp = double(ktmp);
    end
    
    %   Get acquisition size
    [np, nv, ns, nr, ph, zip] = get_key_params(par);
    
    %   Section for reading old USC data
    if strcmp(par.format,'pfile16') && (par.rdb_hdr.user25 > 1)
        
        npts = round(nv*ns / par.rdb_hdr.user25) * par.rdb_hdr.user26;
        ktmp = reshape(ktmp,[np, nv+1, ns, nr]);
        ktmp = permute(ktmp,[1 3 2 4]);
        ktmp = reshape(ktmp,[np, ns*(nv+1), 1, 1, nr]);
        ktmp = ktmp(:,1:npts,:,:,:);
        
    else
        
        %   Reshape data
        if numel(ktmp) == (np*(nv+1)*(zip*ns)*nr)
            rvec1 = [np, nv+1, zip*ns, nr];
            pvec = [1, 3, 2, 4];
            rvec2 = [np, nv*ns, 1, 1, nr];
        elseif numel(ktmp) == (np*(nv+1)*(zip*ns)*nr*ph)
            rvec1 = [np, nv+1, zip*ns, nr, ph];
            pvec = [1, 3, 2, 5, 4];
            rvec2 = [np, nv*ns*ph, 1, 1, nr];
        else
            error('Unexpected data size');
        end
        ktmp = reshape(ktmp,rvec1);
        ktmp = ktmp(:,2:nv+1,:,:,:);    %   Strip out first phase line
        %             ktmp = fft(ktmp,[],3); % This line is only needed if 3DFT is done.
        ktmp = ktmp(:,:,1:ns,:,:);      %   Strip out zip factor
        ktmp = flip(ktmp,1);
        ktmp = permute(ktmp,pvec);
        ktmp = reshape(ktmp,rvec2);
    end
    
    %   Save output
    k = cat(2,k,ktmp); clear ktmp;
    
end
clear ans fid files name path


%   Zero pad fractional echo
if par.image.frac_echo == 1
    np2 = 256;
    k2 = zeros([np2, size(k,2), 1, 1, nr],cls);
    k2(1:np,:,:,:,:) = k;
    k = k2;
    clear k2;
    np = np2;
    par.rdb_hdr.da_xres = np2;
end

%   Read pe table
if exist('GRexttab.txt','file')
    petab = dlmread('GRexttab.txt');
else
    petab = dlmread('GRpetable.txt');
end
petab = petab(1:size(k,2),:);
% petab = petab(6:end,:);

%   Make timestamp
tr = par.image.tr / 1e6;
tstamp = 0:tr:(size(k,2)*tr-tr);


%   Debug. Drop every 2nd point.
% warning('THOWING AWAY DATA!!!');
% k = k(:,1:4:size(k,2),1,1,:);
% petab = petab(1:4:size(petab,1),:);
% tstamp = tstamp(1:4:end);


%   Reduct resolution, if requested
if ~isempty(mtrxL)
    
    %   Define new resolution
    xres = mtrxL(1);
    yres = mtrxL(2);
    zres = mtrxL(3);
    par.rdb_hdr.da_yres = yres+1;
    par.rdb_hdr.da_xres = xres;
    par.image.locsperslab = zres;
    par.rdb_hdr.nslices = par.rdb_hdr.nslices*zres/ns;
    
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
end

%   Reconstruct on larger grid, if requested
if ~isempty(mtrxH)
    xres = mtrxH(1);
    yres = mtrxH(2);
    zres = mtrxH(3);
    par.rdb_hdr.da_yres = yres+1;
    par.rdb_hdr.da_xres = xres;
    par.image.locsperslab = zres;
    par.rdb_hdr.nslices = par.rdb_hdr.nslices*zres/ns;
    
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

    %   Update local variables
    np = xres;
    nv = yres;
    ns = zres;
end

end
