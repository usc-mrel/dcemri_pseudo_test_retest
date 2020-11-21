function [imgR,tform,tformI,Mot] = registervolumes(img,FOV,tform,type,mode)


if nargin < 5 || isempty(mode)
    mode = 'multimodal';
end
if nargin < 4 || isempty(type)
    type = 'rigid';
end
if nargin < 3 || isempty(tform)
    tform = [];
end
if nargin < 2 || isempty(FOV)
    FOV = [1 1 1];
end

%   Open parallel pool
parpl = gcp('nocreate');
if isempty(parpl);
    parpl = parpool(min([size(img,4) feature('numCores')]));
end

%   Matrix size
[~, ~, ~, nt] = size(img);

if isempty(tform)
    %   Estimate motion and apply
    imgR = img;
    img0 = img(:,:,:,end);
    
    if 1
        % parallel
        tform = repmat(affine3d(),[nt,1]);
        parfor i = 1:nt
            [imgR(:,:,:,i),tform(i),Mot(:,i)] = registervolume(img(:,:,:,i),img0,FOV,[],type,mode);
            tformI(i) = invert(tform(i));
        end
    else
        %sequential
        init_tform = affine3d();
        for i=nt:-1:1
            if i < nt
                init_tform = tform(i+1);
            end 
            [imgR(:,:,:,i),tform(i),Mot(:,i)] = registervolume(img(:,:,:,i),img0,FOV,[],type,mode,init_tform);
            tformI(i) = invert(tform(i));
        end
    end
else
    %   Apply only
    imgR = img;
    parfor i = 1:nt
        imgR(:,:,:,i) = registervolume(img(:,:,:,i),[],FOV,tform(i));
        tformI(i) = invert(tform(i));
    end
    Mot = [];
end

end

