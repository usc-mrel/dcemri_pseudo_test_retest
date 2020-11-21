function [kprep,Uprep,nbins] = GRprepkspace(k,pfle,petab,binE,time,type)
%function [kprep,Uprep,nbins] = GRprepkspace(k,pfle,petab,binE,time,type)
%
%   Like GRreshape2 yet combines many preparation steps into the kspace
%   reshape for better memory efficiency
%
%   OBS: the order of k-space samples is slighlty different from how
%   GRreshape does it!
%
%   Yannick 2019

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
[np, nv, ns, nr] = get_key_params2(pfle);
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


kprep = cell(nbins,1);
Uprep = cell(nbins,1);

switch type
    case 0  %   Average all k-space data
        
        parfor b=1:nbins
            kfull = zeros([nv ns np nr],class(k));
            Ufull = zeros([nv ns],class(k));
            for i = find( bin == b )
                y = petab(i,1);
                z = petab(i,2);

                kfull(y,z,:,:) = kfull(y,z,:,:) + reshape(k(i,1,1,:,:), [1 1 np, nr]);
                Ufull(y,z) = Ufull(y,z) + 1;
            end
            
            %   Rearrange kfull
            kfull = permute(kfull,[3 1 2 4]);
            
            %   Normalize data
            Ufull = repmat(reshape(Ufull,[1 nv ns 1]),[np 1 1 1]);
            Uinv = 1./(Ufull+eps);
            Ufull = logical(Ufull);
            for i = 1:nr
                kfull(:,:,:,i) = kfull(:,:,:,i) .* Uinv;
            end
            Ufull = repmat(Ufull,[1 1 1 nr]);
            
            %   Apply fftshift in k-space
            for i = 1:nr
                kfull(:,:,:,i) = fftshiftF(kfull(:,:,:,i),[1 2]);
            end
            
            %   Apply fftshift and convert FE into image domain
            for i = 1:nr
                kfull(:,:,:,i)  = iFastFT(kfull(:,:,:,i),1,1);
                kfull(:,:,:,i)  = fftshift(fftshift(kfull(:,:,:,i),2),3);
                Ufull(:,:,:,i)  = fftshift(fftshift(Ufull(:,:,:,i),2),3);
            end
            
            % store back
            kprep{b} = kfull(Ufull);
            kprep{b} = reshape(kprep{b},[np numel(kprep{b})/np]);
            
            Uprep{b} = [];
            for i = 1:nr
                Uprep{b} = [Uprep{b}(:); find(Ufull(:,:,:,i)) + (b-1)*(nv*ns*np) + (i-1)*(nv*ns*np*nbins)];
            end
        end
    otherwise
        error('Unkown kspace formatting type!')
end
clear k

kprep = cell2mat(reshape(kprep, 1, nbins));
Uprep = cell2mat(reshape(Uprep, nbins, 1));

% sort to make consistent with GRreshape
[Uprep, ind] = sort(Uprep);

kprep = kprep(ind);
kprep = reshape(kprep, [np numel(kprep)/np]);
