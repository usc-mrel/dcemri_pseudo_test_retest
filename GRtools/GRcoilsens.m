function [S,sf] = GRcoilsens(k,NCAL)
%   Function to compute+save or load coil sensitivity

if nargin < 2 || isempty(NCAL)
    NCAL = 40^2;
end

%%  Get parameters
np = size(k,1);
nv = size(k,2);
ns = size(k,3);
nr = size(k,5);


if (1)
    %%  Estimate sensitivity
    fr = 0.9*sqrt(NCAL);
    fw = 0.015;
    F = FermiFilt([np nv ns],[2*fr fr fr],fw);
    img = iFastFT(k(:,:,:,1,:) .* repmat(F,[1 1 1 1 nr]),[1 2 3],1);
    sf = 1/max(abs(img(:)));
    S = img ./ repmat(sqrt(sum(abs(img).^2,5)),[1 1 1 1 nr]);
    
else
    
    %%  eSPRIRIT
    img = iFastFT(k,[1 2 3],1);
    sf = 1/max(abs(img(:)));
    clear img;
    kcal = extract_kcal(k(:,:,:,1,:),[sqrt(NCAL) sqrt(NCAL) sqrt(NCAL)]);
    opt = SPSENSE_optset;
    opt.kernel = [3 7 7];
    opt.kernreg = 'auto';
    opt.size = [sqrt(NCAL) sqrt(NCAL) sqrt(NCAL) 1 nr];
    G = gKERN3d(kcal,opt);
    S = sMAP(G,opt);
    for i = 1:nr
        S2(:,:,:,1,i) = imresize3d(S(:,:,:,1,i),[np nv ns]);
    end
    S = S2;clear S2;
    clear G opt kcal
end

end




%% ACTUAL SUB FUNCTIONS
function H = FermiFilt(sz,N,T)

%   Dim 1
E = N(1)/sz(1);
x = -1:2/(sz(1)-1):1;
H1 = 1./(1+exp((abs(x)-E)/T))';
H = repmat(H1,[1,sz(2),sz(3)]);

%   Dim 2
E = N(2)/sz(2);
x = -1:2/(sz(2)-1):1;
H1 = reshape(1./(1+exp((abs(x)-E)/T)),[1 sz(2) 1]);
H = H.*repmat(H1,[sz(1),1,sz(3)]);

%   Dim3
E = N(3)/sz(3);
x = -1:2/(sz(3)-1):1;
H1 = reshape(1./(1+exp((abs(x)-E)/T)),[1 1 sz(3)]);
H = H.*repmat(H1,[sz(1),sz(2),1]);

end