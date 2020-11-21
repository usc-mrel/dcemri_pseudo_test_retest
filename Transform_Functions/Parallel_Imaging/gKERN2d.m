function G = gKERN2d(kcal,opt)
%   Computes the 2D GRAPPA/SPIRiT kernel
%   
%   Author: RML
%   Date: 02/2011
%   
%   Usage: G = gKERN2d(kcal,opt)
%   
%   Input:
%   kcal: Calibration k-space of size RO x PE x NS x NT x NR
%   opt: Specifies transform options for recon. Use SPSENSE_optset.m
%   
%   Output:
%   G: GRAPPA/SPIRiT kernel


%   Get sizes
[~, ~, ns , nt, nr] = size(kcal);
kf = opt.kernel;

%   Loop through slices (outer loop since totally independent)
G = zeros([kf ns nt nr nr],opt.class);
for ins = 1:ns
    
    %   Loop through reciver coils, construct the neighbourhood matrix (N)
    N = [];
    for inr = 1:nr
        N = [N im2col(squeeze(kcal(:,:,ins,1,inr)),opt.kernel,'sliding').'];
    end
    
    %   Solve for GRAPPA weights for each coil
    for inr = 1:nr
        %   Create kernel sampling array and indices
        samp = ones([kf nr],opt.class);
        samp((end+1)/2,(end+1)/2,inr) = 0;  %   Exclude points to be estimated
        curr = ~samp;                       %   Include points to be estimated
        indS = find(samp);                  %   Define sampled indices
        indC = find(curr,1);                %   Define current location index
        
        %   Remove the current location from the neighbourhood matrix
        %   Create the position vector
        X = N(:,indC);
        NS = N(:,indS);
        
        %   Define some intermediate variables
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
        
        %   Store GRAPPA kernel
        G(:,:,ins,1,:,inr) = flipdim(flipdim(Gtmp,1),2);
    end
    
    %   Display progress
    if opt.verbose
        fprintf(' %u%%',round(ins/ns*100));
    end
end
clear Gtmp NS N NSt X curr samp indC indS ins inr lambda

%   Pad G and convert to image space
%   (multiply in image space rather than convolve in k-space)
Gt = G * opt.size(1)*opt.size(2);
G = zeros([opt.size(1) opt.size(2) opt.size(3) 1 nr nr],opt.class);
G(ceil((opt.size(1)-kf(1))/2)+1:ceil((opt.size(1)+kf(1))/2),...
    ceil((opt.size(2)-kf(2))/2)+1:ceil((opt.size(2)+kf(2))/2),:,:,:,:) = Gt;
clear Gt
opt.FTshift = 1;
opt.FTdim = [1 2];
G = opt.iFT(G,opt);

end

%   Function to estimate regularization factor
function reg = est_reg(NSt,NS,X,opt)

    %   Define array of possible regularizaiton factors
    N = 12;
    mn = 1e-6;
    mx = 1e6;
    reg = exp(log(mn):(log(mx)-log(mn))/N:log(mx));
    reg2 = exp(log(mn/100):(log(mx*100)-log(mn/100))/(100*N):log(mx*100));
    
    %   Loop through to calculate fitting error
    scale = norm(NSt,'fro')/size(NSt,1);
    err = zeros([1,N+1],opt.class);
    for indx = 1:N+1
        lambda = scale * reg(indx);
        GR = pinv(NSt +  lambda*eye(size(NSt)))*NS' * X;
        dif = NS*GR - X;
        err(indx) = norm(dif(:));
    end
    clear GR dif scale indx lambda NS NSt X
    
    %   Rescale for fitting
%     hold on;
%     plot(log(reg),err,'x-');
    lr = log(reg);
    lr2 = log(reg2);
%     le = log(err);
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
