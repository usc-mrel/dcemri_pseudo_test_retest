function [img,opt] = SP_STATIC_MOCO(k,opt)

%%  Setup
opt.worder = calc_wave_order(ones(opt.size(1:3)),16);
nt = opt.size(4);
opt.phase_corr = [1 2 3];
DC = [];

%   Create initial undersampled image
%fprintf(opt.fid,'\tEstimating initial image...');tic
img = iFFTU(k,opt);
img = iSENSE(img,opt);
img = mean(img,4);
img = repmat(img,[1 1 1 nt 1]);
%fprintf(opt.fid,'done (%g sec)\n',toc);

%%  Initial sparse sense to get a half decent image
for i = 1:5
    
    %%  Replace data
    %   Convert to k-space
    img = fSENSE(img,opt);
    img = fFastFT(img,opt);
    
    %   Replace known data
    DC = [DC sqrt(sum(abs(img(opt.U) - k).^2))];
    img(opt.U) = k;
    
    %   Convert to coil-combined dynamic image
    img = iFastFT(img,opt);
    img = iSENSE(img,opt);
    
    %%  Dynamic combine
    %   Convert to static image
    img = mean(img,4);
    
    %   Wavelet
    WC = fcwtN(img,'',opt.worder);
    [WC, lambda] = CELLsoftthresh(WC,opt.lambda,1);
    img = icwtN(WC,'',opt.worder);clear WC;
    
    %%  Plotting
    figure(1);
    semilogy(1:length(DC),DC./DC(1),'x-');
    grid on;
    drawnow;
    figure(2);
    images(abs(img(:,:,60)));
    drawnow;
    
    %   Dynamic expansion
    img = repmat(img,[1 1 1 nt+1 1]);
    
end

%%  POCS iterations
tform = [];
tformI = [];
%fprintf(opt.fid,'\tCore POCS iterations:\n');tic
mxit = 150;
for i = 0:mxit
    i
    
    %%  Replace data
    %   Convert to k-space
    img = fSENSE(img,opt);
    img = fFastFT(img,opt);
    
    %   Replace known data
    DC = [DC sqrt(sum(abs(img(opt.U) - k).^2))];
    img(opt.U) = k;
    
    %   Convert to coil-combined dynamic image
    img = iFastFT(img,opt);
    img = iSENSE(img,opt);
    
    
    %%  Motion correction
    %   Remove phase
    opt.Ph = exp(sqrt(-1)*angle(img));
    img = real(opt.iPH(img,opt));
    
    %   Estimate/apply motion correction
    img = cat(4,mean(img,4),img);
    if mod(i,5) == 0
        [img,tform,tformI,Mot] = registervolumes(img,opt.FOV);
        figure(4);
        plot(1:size(img,4),Mot);
        drawnow;
    elseif i>=5
        img = registervolumes(img,opt.FOV,tform);
    end
    img = img(:,:,:,2:end);
    
    %%  Dynamic combine
    %   Convert to static image
    img = mean(img,4);
    opt.Ph = mean(opt.Ph,4);
    
    %%   Image sparsity
    WC = fcwtN(img,'',opt.worder);
    [WC, lambda] = CELLsoftthresh(WC,opt.lambda,1);
    img = icwtN(WC,'',opt.worder);clear WC;
    img = real(img);
    
    %%   Terminate
    if i == mxit
        img = opt.fPH(img,opt);
        return;
    end
    
    
    %%  Plotting
    figure(1);
    semilogy(1:length(DC),DC./DC(1),'x-');
    grid on;
    drawnow;
    figure(2);
    images(img(:,:,60));
    drawnow;
    

    %%   Un-register
    img = opt.fPH(img,opt);
    img = repmat(img,[1 1 1 nt+1 1]);
    img = registervolumes(img,opt.FOV,tformI);
    img = img(:,:,:,2:end);

    
    %fprintf(opt.fid,'\t\t%g: lambda = %f\n',i,lambda);

end
%fprintf(opt.fid,'done (%g sec)\n',toc);

img;


