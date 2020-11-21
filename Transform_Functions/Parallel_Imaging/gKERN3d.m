function G = gKERN3d(kcal,opt)
%   Computes the 3D GRAPPA/SPIRiT kernel
%   
%   Author: RML
%   Date: 02/2011
%   
%   Usage: G = gKERN3d(kcal,opt)
%   
%   Input:
%   kcal: Calibration k-space of size RO x PE x PE2 x NT x NR
%   opt: Specifies transform options for recon. Use SPSENSE_optset.m
%   
%   Output:
%   G: GRAPPA/SPIRiT kernel


%   Get sizes
[ro pe ns nt nr] = size(kcal);
kf = opt.kernel;
indro = ceil(kf(1)/2):ro-floor(kf(1)/2);
indpe = ceil(kf(2)/2):pe-floor(kf(2)/2);
indns = ceil(kf(3)/2):ns-floor(kf(3)/2);
hkf = floor(kf/2);
%   Loop through reciver coils, construct the neighbourhood matrix (N)
N = [];
for inr = 1:nr
    Nt = zeros([length(indro)*length(indpe)*length(indns),prod(kf)],opt.class);
    ind = 1;
    for iro = indro
    for ipe = indpe
    for ins = indns
        temp = kcal(iro-hkf(1):iro+hkf(1),ipe-hkf(2):ipe+hkf(2),ins-hkf(3):ins+hkf(3),1,inr);
        Nt(ind,:) = temp(:).';
        ind = ind+1;
    end
    end
    end
    N = [N Nt];
end
clear Nt indro indpe indns ind iro ipe ins temp

%   Solve for GRAPPA weights for each coil
G = zeros([kf nt nr nr],opt.class);
for inr = 1:nr
    %   Create kernel sampling array and indices
    samp = ones([kf nr],opt.class);
    samp((end+1)/2,(end+1)/2,(end+1)/2,inr) = 0;    %   Exclude points to be estimated
    curr = ~samp;                                   %   Include points to be estimated
    indS = find(samp);
    indC = find(curr,1);
    
    %   Remove the current location from the neighbourhood matrix
    %   Create the position vector
    X = N(:,indC);
    NS = N(:,indS);
    Gtmp = zeros([kf nr],opt.class);
    NSt = NS'*NS;
    
    %   Estimate regularization factor
    %   ---------EXPERIMENTAL---------
    if strcmp(opt.kernreg,'auto');
        opt.kernreg = est_reg(NSt,NS,X,opt);
    end
    lambda = norm(NSt,'fro')/size(NSt,1) * opt.kernreg;
    
    %   Invert to solve for GRAPPA weights
    Gtmp(indS) = pinv(NSt +  lambda*eye(size(NSt)))*NS' * X;
%     Gtmp(indS) = (NSt +  lambda*eye(size(NSt)))\NS' * X;
    
    %   Store GRAPPA kernel
    G(:,:,:,1,:,inr) = flip(flip(flip(Gtmp,1),2),3);
    
    %   Display progress
    if opt.verbose
        fprintf(' %g',inr);
    end
end

%   Define downsampling
%downs = [2 2 1];
downs = [1 1 1];
Gs = ceil([opt.size(1)/downs(1) opt.size(2)/downs(2) opt.size(3)/downs(3)]);

%   Pad G and convert to image space
%   (multiply in image space rather than convolve in k-space)
if opt.kernel(1) == 1
    G = G * Gs(2)*Gs(3);
    Gf = zeros([1 Gs(2) Gs(3) 1 nr nr],opt.class);
    Gf(1,ceil((Gs(2)-kf(2))/2)+1:ceil((Gs(2)+kf(2))/2),...
        ceil((Gs(3)-kf(3))/2)+1:ceil((Gs(3)+kf(3))/2),:,:,:) = G;
    opt.FTshift = 1;
    G = opt.iFT(Gf,opt);
    G = repmat(G,[Gs(1) 1 1 1 1]);
else
    G = G * Gs(1)*Gs(2)*Gs(3);
    Gf = zeros([Gs(1) Gs(2) Gs(3) 1 nr nr],opt.class);
    Gf(ceil((Gs(1)-kf(1))/2)+1:ceil((Gs(1)+kf(1))/2),...
       ceil((Gs(2)-kf(2))/2)+1:ceil((Gs(2)+kf(2))/2),...
       ceil((Gs(3)-kf(3))/2)+1:ceil((Gs(3)+kf(3))/2),:,:,:) = G;
    opt.FTshift = 1;
    opt.FTdim = [1 2 3];
    for i = 1:nr
        for j = 1:nr
            Gf(:,:,:,:,i,j) = opt.iFT(Gf(:,:,:,:,i,j),opt);
        end
    end
    G = Gf;
    clear Gf;
end

end


%   Function to estimate regularization factor
function reg = est_reg(NSt,NS,X,opt)

    %   Define array of possible regularizaiton factors
    N = 12;
    mn = 1e-4;
    mx = 1e4;
    reg = exp(log(mn):(log(mx)-log(mn))/N:log(mx));
    reg2 = exp(log(mn/100):(log(mx*100)-log(mn/100))/(100*N):log(mx*100));
    
    %   Loop through to calculate fitting error
    scale = norm(NSt,'fro')/size(NSt,1);
    err = zeros([1,N+1],opt.class);
    for indx = 1:N+1
        lambda = scale * reg(indx);
        GR = pinv(NSt +  lambda*eye(size(NSt)))*NS' * X;
%         GR = (NSt +  lambda*eye(size(NSt)))\NS' * X;
        dif = NS*GR - X;
        err(indx) = norm(dif(:));
    end
    clear GR dif scale indx lambda NS NSt X
    
    %   Rescale for fitting
    lr = log(reg);
    lr2 = log(reg2);
    le = err;
    le = le-min(le);
    le = le./max(le);
    
    %   Fit with sigmoid function
    b = glmfit(lr,[le' ones([N+1,1],opt.class)],'binomial','link','probit');
    y = glmval(b,lr2,'probit','size',ones([length(lr2),1],opt.class));
    if opt.plot
        semilogx(reg,le,'x',reg2,y);axis tight;grid on;
        ylabel('Kernel Error');xlabel('Regularization Factor');drawnow;
    end
    
    %   Find theshold value
    %   This assumes that the initial deflection of the sigmoid represents
    %   the optimal regularization parameter. Untested.
    th = 0.01;
    reg = reg2(find(abs(y-th) == min(abs(y-th)),1));
    
    %   Display value
    if opt.verbose
        fprintf('Regularization set to %g...',reg);
    end
    
    
end
