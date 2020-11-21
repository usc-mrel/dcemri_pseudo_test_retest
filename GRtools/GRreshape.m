function [kfull,Ufull,kwt] = GRreshape(k,petab,par,binE,time,type)

if nargin < 6 || isempty(type)
    type = 0;
end
if nargin < 4 || isempty(binE)
    binE = Inf;
end
if nargin < 3
    error('need at least 3 inputs');
end

%   Get acquisition parameters
[np, nv, ns, nr] = get_key_params(par);
nviews = size(k,2);

%   Define bin edges
if isinf(binE)
    binE = time(end);
end
if binE(end) < time(end)
    binE(end) = time(end);
end

%   Put samples into bins
bin = zeros(1,nviews);
for i = 1:nviews
    bin(i) = find(binE-time(i)>=0,1);
end
nbins = max(bin);

%   Rearrange k for better memory access
k = permute(k,[2 3 4 1 5]);

kfull = zeros([nv ns nbins np nr],class(k));
Ufull = zeros([nv ns nbins],class(k));
switch type
    case 0  %   Average all k-space data
        for i = 1:nviews
            y = petab(i,1);
            z = petab(i,2);
            t = bin(i);
            
            kfull(y,z,t,:,:) = kfull(y,z,t,:,:) + k(i,1,1,:,:);
            Ufull(y,z,t) = Ufull(y,z,t) + 1;
        end
    case 1  %   Use early points
        for i = nviews:-1:1
            y = petab(i,1);
            z = petab(i,2);
            t = bin(i);
            
            kfull(y,z,t,:,:) = k(i,1,1,:,:);
            Ufull(y,z,t) = 1;
        end
        
    case 2  %   Use late points
        for i = 1:nviews
            y = petab(i,1);
            z = petab(i,2);
            t = bin(i);
            
            kfull(y,z,t,:,:) = k(i,1,1,:,:);
            Ufull(y,z,t) = 1;
        end
    case 3  %   Only for three frames: [late; average; early]
        if nbins == 3
            %   Frame 1 (favour late views)
            for i = find(bin==1,1,'first'):find(bin==1,1,'last')
                y = petab(i,1);
                z = petab(i,2);
                t = bin(i);
                
                kfull(y,z,t,:,:) = k(i,1,1,:,:);
                Ufull(y,z,t) = 1;
            end
            %   Frame 2 (Average)
            for i = find(bin==2,1,'first'):find(bin==2,1,'last')
                y = petab(i,1);
                z = petab(i,2);
                t = bin(i);
                
                kfull(y,z,t,:,:) = kfull(y,z,t,:,:)+k(i,1,1,:,:);
                Ufull(y,z,t) = Ufull(y,z,t)+1;
            end
            %   Frame 3 (favour early views)
            for i = find(bin==3,1,'last'):-1:find(bin==3,1,'first')
                y = petab(i,1);
                z = petab(i,2);
                t = bin(i);
                
                kfull(y,z,t,:,:) = k(i,1,1,:,:);
                Ufull(y,z,t) = 1;
            end
        end
end
clear k

%   Rearrange kfull
kfull = permute(kfull,[4 1 2 3 5]);

%   Normalize data
Ufull = repmat(reshape(Ufull,[1 nv ns nbins 1]),[np 1 1 1 1]);
Uinv = 1./(Ufull+eps);
if nargout == 3
    kwt = sqrt(Ufull);
end
Ufull = logical(Ufull);
for i = 1:nr
    kfull(:,:,:,:,i) = kfull(:,:,:,:,i) .* Uinv;
end
Ufull = repmat(Ufull,[1 1 1 1 nr]);

if nargout == 3
    kwt = repmat(kwt,[1 1 1 1 nr]);
end

%   Apply fftshift in k-space
for i = 1:nr
    kfull(:,:,:,:,i) = fftshiftF(kfull(:,:,:,:,i),[1 2]);
end
