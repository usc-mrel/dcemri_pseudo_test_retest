function GRwrite_dicom(img,pfile,header,subdir,namebase,se_desc,se_off,install)

%   Check input to export dicoms back to scanner
if nargin<6
    install = 0;
end

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


if nargin < 7
    se_off = 0;
end
if nargin < 6
    se_desc = header.SeriesData.se_desc;
end
if nargin < 5
    namebase = 'Image';
end

%   Check that image is single
if ~isa(img,'single')
    img = single(img);
end

%   Get some image info
[np,nv,ns,nt] = size(img);
xr = header.RawHeader.rc_xres;
yr = header.RawHeader.rc_yres;
zr = header.RawHeader.nslices/header.RawHeader.npasses;

%   Check that image matches pfile
if (xr ~= np) || (yr ~= nv) || (zr ~= ns)
    img2 = zeros([xr, yr, zr, nt],'single');
    parfor i = 1:nt
    %for i = 1:nt
    
        %   Shrink dimensions
        imgtmp = zeros([min([np nv ns],[xr yr zr]) nt],'like',img);
        if any([np nv ns] > [xr yr zr])
            imgtmp(:,:,:,i) = real(imshrinkFT(img(:,:,:,i),min([np nv ns],[xr yr zr]),Inf));
        else
            imgtmp(:,:,:,i) = img(:,:,:,i);
        end
        
        %   Interpolate dimensions
        if any([np nv ns] < [xr yr zr])
            img2(:,:,:,i) = real(imresizeFT(imgtmp(:,:,:,i),[xr, yr, zr],header));
        else
            img2(:,:,:,i) = imgtmp(:,:,:,i);
        end
    end
    img = img2;
    [np,nv,ns,nt] = size(img);
    clear img2 imgtemp
end


%   Get corner points and orientation
for s = 1:ns
    % Get corners and orientation for this slice location
    corners(s) = GERecon('Pfile.Corners', s); %#ok<AGROW>
    orientation(s) = GERecon('Pfile.Orientation', s); %#ok<AGROW>
end


%   Rotate/transpose
trans = header.RawHeader.transpose;
rot   = header.RawHeader.rotation;
parfor p = 1:nt
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


%   Gradwarp
parfor p = 1:nt
%for p = 1:nt
    for s = 1:ns
        img(:,:,s,p) = GERecon('Gradwarp',img(:,:,s,p),corners(s));
    end
end


%   Remove negative values
%img(img<0) = 0;

%   Make a few custom tags
SeNum = 100*header.SeriesData.se_no + se_off;
% customTag.Group = hex2dec('0028');
% customTag.Element = hex2dec('0030');
% customTag.VRType = 'DS';
% customTag.Value = [num2str(header.RawHeader.fov/np) ' / ' num2str(header.RawHeader.fov/nv)];

wdwW.Group = hex2dec('0028');
wdwW.Element = hex2dec('1051');
wdwW.VRType = 'DS';
wdwL.Group = hex2dec('0028');
wdwL.Element = hex2dec('1050');
wdwL.VRType = 'DS';
for p = 1:nt
for s = 1:ns
    
    % Orient the image
    finalImage = GERecon('Orient', img(:,:,s,p), orientation(s));
    finalCorners = GERecon('Orient', corners(s), orientation(s));
    
    %   Negative pixels don't look nice
    wdwW.Value = num2str(round(max(finalImage(:))));
    wdwL.Value = num2str(round(max(finalImage(:))/2));
    
    % Save DICOMs
    imageNumber = ImageNumber(s, 1, p, pfile, header);
    filename = ['Ox_DICOM/' subdir '/' namebase '_P' num2str(p-1,'%03u') '_S' num2str(s-1,'%03u') '.dcm'];
    GERecon('Dicom.Write', filename, finalImage, imageNumber, orientation, finalCorners,...
            SeNum,se_desc,wdwW,wdwL);
end
end



%   install images to scanner database
if install
    
    hid = fopen('hostname.txt');
    if hid
        hostnme = textscan(hid,'%s');
        fclose(hid);
        hostnme = hostnme{1};
        hostnme = hostnme{1};
        cmdbase = ['scp -r Ox_DICOM/' subdir];
        switch hostnme
            case 'VMopenSUSE'
                host = 'mlebel@172.22.37.180';
                rdir = ['/home/mlebel/Desktop/Ox_' subdir];
                cmd = [cmdbase ' ' host ':' rdir];
                system(cmd);
                cmd2 = ['ssh ' host ' "mv ' rdir ' ' rdir '.sdc_open"'];
                system(cmd2);
                
            case 'MRGECAA1'
                host = 'sdc@172.22.13.7';
                rdir = ['/export/home1/sdc_image_pool/import/Ox_' subdir];
                cmd = [cmdbase ' ' host ':' rdir];
                system(cmd);
                cmd2 = ['ssh ' host ' "mv ' rdir ' ' rdir '.sdc_open"'];
                system(cmd2);
                
            case 'MRGESEA1'
                host = 'sdc@139.48.44.90';
                rdir = ['/export/home1/sdc_image_pool/import/Ox_' subdir];
                cmd = [cmdbase ' ' host ':' rdir];
                system(cmd);
                cmd2 = ['ssh ' host ' "mv ' rdir ' ' rdir '.sdc_open"'];
                system(cmd2);
                
            otherwise
                warning('Unknown host');
        end
        
    else
        warning('Unable to check source host');
    end
    
end

end

function number = ImageNumber(slice, echo, phase, pfile, header)
% Image numbering scheme:
% P0S0E0, P0S0E1, ... P0S0En, P0S1E0, P0S1E1, ... P0S1En, ... P0SnEn, ...
% P1S0E0, P1S0E1, ... PnSnEn
    slicesPerPhase = header.RawHeader.nslices * pfile.echoes;
    number = (phase-1) * slicesPerPhase + (slice-1) * pfile.echoes + (echo-1);
end
