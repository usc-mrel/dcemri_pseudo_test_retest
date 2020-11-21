function [R,shifts] = phasematch3D(k,dim)

%   Check inputs
if nargin < 2
    dim = '123';
end

%   Define some variables
us = 1/100;
interp = 'spline';

%   Get size
[np nv ns nt nr] = size(k);

%   Create output array
R = ones([np nv ns nt nr]);
shifts  = zeros(nt,3);

%   Loop through time dimensions
for int=2:nt
    
    %   Compute cross correlation
    C = k(:,:,:,1,:).*conj(k(:,:,:,int,:));
    C = C./(abs(C)+eps);
    C = abs(iFFT(C,[1 2 3],0));
    C = sum(C,5);
    C = fftshift(C);
    
    %   Locate peak
    ind = find(C(:) == max(C(:)),1);
    [i,j,l] = ind2sub([np nv ns],ind);
        
    %   Update phase shift arrary for vertical then horizontal then through plane
    if ~isempty(strfind(dim,'1'))
        %   Find first direction shift
        r = interp1(1:np,C(:,j,l),us:us:np,interp);
        shifts(int,1) = mean(find(r==max(r))*us) - ceil((0.5+np)/2);
        
        ph_ramp = exp(-sqrt(-1)*2*pi*shifts(int,1)/2*(-1:2/np:1-1/np))';
        R(:,:,:,int,:) = R(:,:,:,int,:).*repmat(ph_ramp,[1 nv ns 1 nr]);
    end
    if ~isempty(strfind(dim,'2'))
        %   Find second direction shift
        r = interp1(1:nv,C(i,:,l),us:us:nv,interp);
        shifts(int,2) = mean(find(r==max(r))*us) - ceil((0.5+nv)/2);
        
        ph_ramp = exp(-sqrt(-1)*2*pi*shifts(int,2)/2*(-1:2/nv:1-1/nv));
        R(:,:,:,int,:) = R(:,:,:,int,:).*repmat(ph_ramp,[np 1 ns 1 nr]);
    end
    if ~isempty(strfind(dim,'3'))
        %   Find second direction shift
        r = interp1(1:ns,squeeze(C(i,j,:)),us:us:ns,interp);
        shifts(int,3) = mean(find(r==max(r))*us) - ceil((0.5+ns)/2);
        
        ph_ramp = exp(-sqrt(-1)*2*pi*shifts(int,2)/2*(-1:2/ns:1-1/ns));
        R(:,:,:,int,:) = R(:,:,:,int,:).*repmat(reshape(ph_ramp,[1 1 ns]),[np nv 1 1 nr]);
    end
end
