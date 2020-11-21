function [par, img] = dicomr(path,dt)
%   Looping function to read in multiple dicom files
%   Required functions: dicomread, dicominfo (image proc. toolbox)
%   
%   Usage: [par, img] = dicomr(path)
%   Author: R. Marc Lebel
%   Date: 04/05/2006
%   
%   Input:
%   fin: (optional) path to a directory containing dicom (.IMA) files.
%   dt: flag for data typeL 'd' double, 's' single, 'i' uint16
%   
%   Output:
%   par: the dicom meta data from the first file
%   img: stack of images (in double format)

%   Use current directory as default
if nargin<1 || isempty(path)
    path = cd;
end
cd(path);

if nargin < 2
    dt = 'd';
end

%   Obtain a directory listing
files = dir;
nfiles = size(files);
nfiles = nfiles(1);

%   Loop through files, determine if they are dicom, then read them in
count = 1;
for i = 1:nfiles
    name = files(i).name;
    if length(name) > 4
        if ~strcmpi(name(1),'.') && ...
           ~strcmpi(name(end-3:end),'.txt') && ...
           ~strcmpi(name(end-3:end),'.nii') && ...
           (strcmpi(name(end-3:end),'.IMA') || ...
            strcmpi(name(end-3:end),'.dcm') || ...
            ~isempty(regexp(name,'MRDC','match')) || ...
            ~isempty(regexp(name,'IM-','match')))
%             contains(name,'MRDC') || ... % guess what this does not work on clipper
%             contains(name,'IM-'))
            
            %   Obtain parameters
            partemp = dicominfo(name);
            
            %   Get parameters from first image
            if count == 1
                par = partemp;
            end
            
            %   Update specific parameters
            par.SliceLocation(:,count) = partemp.SliceLocation(:);
            par.ImageOrientationPatient(:,count) = partemp.ImageOrientationPatient(:);
            par.ImagePositionPatient(:,count) = partemp.ImagePositionPatient(:);
            par.InstanceNumber(:,count) = partemp.InstanceNumber(:);
            par.EchoTime(:,count) = partemp.EchoTime(:);
            
            %   Read in image
            if strcmp(dt,'i')
                img(:,:,count) = uint16(dicomread(name));
            elseif strcmp(dt,'s')
                img(:,:,count) = single(dicomread(name));
            else
                img(:,:,count) = double(dicomread(name));
            end
            
            %   Increment array size counter
            count=count+1;
        end
    end
end

%   Sort Images
[~,ind] = sort(par.InstanceNumber);
% [~,ind] = sort(par.SliceLocation);
img = img(:,:,ind);
par.SliceLocation = par.SliceLocation(ind);
par.ImagePositionPatient = par.ImagePositionPatient(:,ind);
par.ImageOrientationPatient = par.ImageOrientationPatient(:,ind);
par.InstanceNumber = par.InstanceNumber(ind);

return
