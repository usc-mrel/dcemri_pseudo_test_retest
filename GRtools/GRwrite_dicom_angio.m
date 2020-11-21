function GRwrite_dicom_angio(img,namebase,se_desc,se_off)

%   Define subdir
subdir = 'Angio';

%   Check that directory exists and is empty
if exist('Ox_DICOM','dir')
    if exist(['Ox_DICOM/' subdir],'dir')
        delete(['Ox_DICOM/' subdir '/*']);
    else
        mkdir(['Ox_DICOM/' subdir]);
    end
else
    mkdir(['Ox_DICOM/' subdir]);
end


%   Load Pfile info
pfiles = dir([pwd '/P*.7']);
fname = pfiles(1).name;
%pfile = GERecon('Pfile.Load',fname,'No-Anonymize');
pfile = GERecon('Pfile.Load',fname);
header = GERecon('Pfile.Header',pfile);


if nargin < 4
    se_off = 0;
end
if nargin < 3
    se_desc = header.SeriesData.se_desc;
end
if nargin < 2
    namebase = 'Angio';
end

%   Check that image is single
if ~isa(img,'single')
    img = single(img);
end

%   Get some image info
[np,nv,ns,nt] = size(img);
xr = header.RawHeader.rc_xres;
yr = header.RawHeader.rc_yres;
zr = pfile.slices;

%   Check that image matches pfile
if (xr ~= np) || (yr ~= nv)
    error('Wrong image size');
end



%   Get corner points and orientation for center
corners = GERecon('Pfile.Corners', round(zr/2));
orientation = GERecon('Pfile.Orientation', round(zr/2));


%   Rotate/transpose
trans = header.RawHeader.transpose;
rot   = header.RawHeader.rotation;
for p = 1:nt
for s = 1:ns
    if trans > 0
        img(:,:,s,p) = img(:,:,s,p)';
    end
    if rot > 0
        img(:,:,s,p) = rot90(img(:,:,s,p),rot);
    end
    if trans < 0
        img(:,:,s,p) = img(:,:,s,p)';
    end
end
end


%   Remove negative values
img(img<0) = 0;

%   Make a few custom tags
SeNum = 100*header.SeriesData.se_no + se_off;


for p = 1:nt
for s = 1:ns
    
    % Orient the image
    finalImage = GERecon('Orient', img(:,:,s,p), orientation);
    finalCorners = GERecon('Orient', corners, orientation);
    
    % Save DICOMs
    imageNumber = ImageNumber(s, 1, p, pfile);
    filename = ['Ox_DICOM/' subdir '/' namebase '_P' num2str(p-1,'%03u') '_S' num2str(s-1,'%03u') '.dcm'];
    GERecon('Dicom.Write', filename, finalImage, imageNumber, orientation, finalCorners,SeNum,se_desc);
end
end


end

function number = ImageNumber(slice, echo, phase, pfile)
% Image numbering scheme:
% P0S0E0, P0S0E1, ... P0S0En, P0S1E0, P0S1E1, ... P0S1En, ... P0SnEn, ...
% P1S0E0, P1S0E1, ... PnSnEn
    slicesPerPhase = pfile.slices * pfile.echoes;
    number = (phase-1) * slicesPerPhase + (slice-1) * pfile.echoes + (echo-1);
end
