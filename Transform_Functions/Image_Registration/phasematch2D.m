function [R,shifts] = phasematch2D(k,dim)

%   Check inputs
if nargin < 2
    dim = '12';
end

%   Define some variables
us = 1/100;
interp = 'spline';

%   Get size
[np nv ns nt nr] = size(k);

%   Create output array
R = ones([np nv ns nt nr]);
shifts  = zeros(nt,2);

%   Loop through outter dimensions
for int=2:nt
    
    %   Compute cross correlation (averaged accross slices and receivers)
    C = k(:,:,:,1,:).*conj(k(:,:,:,int,:));
    C = C./(abs(C)+eps);
    C = abs(ift2d(C,1));
    C = sum(sum(C,3),5);
    
    %   Locate peak
    ind = find(C(:) == max(C(:)),1);
    [i,j] = ind2sub([np nv],ind);
    
    %   Update phase shift arrary for vertical then horizontal
    if ~isempty(strfind(dim,'1'))
        %   Find first direction shift
        r = interp1(1:np,C(:,j),us:us:np,interp);
        shifts(int,1) = mean(find(r==max(r))*us) - ceil((0.5+np)/2);
        
        %   Make phase ramp
        ph_ramp = exp(-sqrt(-1)*2*pi*shifts(int,1)/2*(-1:2/np:1-1/np))';
        R(:,:,:,int,:) = R(:,:,:,int,:).*repmat(ph_ramp,[1 nv ns 1 nr]);
    end
    
    if ~isempty(strfind(dim,'2'))
        %   Find second direction shift
        r = interp1(1:nv,C(i,:),us:us:nv,interp);
        shifts(int,2) = mean(find(r==max(r))*us) - ceil((0.5+nv)/2);
        
        %   Make phase ramp
        ph_ramp = exp(-sqrt(-1)*2*pi*shifts(int,2)/2*(-1:2/nv:1-1/nv));
        R(:,:,:,int,:) = R(:,:,:,int,:).*repmat(ph_ramp,[np 1 ns 1 nr]);
    end
end
