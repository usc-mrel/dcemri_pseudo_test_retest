function tbolus = GRestbolus(k0,tstamp0)
%   Locates the arrival of contast agent using the center of k-space
%
%   Authour: R. Marc Lebel
%   Date: March 9 2016
%   
%   Inputs (optional):
%   k0: x-t space data with N readout points (image space) x T time points. Time points do not need
%       to be evenly spaced
%   tstamp0: time stamps corresponding to the temporal locations of the k0 samples
%   
%   Output:
%   tbolus: Estimated arrival time of the contrast bolus.

%   Debug/vis flag (useful for seeing fits and for exporting data for plotting elsewhere)
export = false;
plot_fit = true;

%   Read data if not provided.
if nargin ~= 2
    %   Get centre of k-space
    %   NEED TO KNOW NUMBER OF READOUT POINTS TO REQUEST
    [k0,~,~,tstamp0] = GRread2('./',[0 1 1]);
    %[k0,~,~,tstamp0] = GRread2('./',[1 1 1]);
    k0 = ifftshift(ifft(ifftshift(k0,1),[],1),1);
    k0 = sqrt(sum(abs(k0).^2,5));
end
k0 = double(k0(:,:,:));


%   Filter to remove spikes (likely due to lack of disdaqs during multi-phase scan)
ks = k0;
for i = 1:size(k0,1);
    ks(i,:) = medfilt1(k0(i,:),3);
end

%   Interpolate k to uniform spacing
dt = 1.0; % seconds
tstamp = min(tstamp0):dt:max(tstamp0);
k = zeros(size(k0,1),length(tstamp));
for i = 1:size(k0,1);
    k(i,:) = interp1(tstamp0,ks(i,:),tstamp);
end
tstamp = tstamp(:);
clear ks;

%   Perform linear filtering then estimate the goodness of fit
kf = k;
R2 = zeros(size(k,1),1);
for i = 1:size(k,1)
    kf(i,:) = smooth(k(i,:),9);
    SSres = sum(abs(k(i,:) - kf(i,:)).^2);
    SStot = sum(abs(k(i,:) - mean(k(i,:))).^2);
    R2(i) = 1-SSres/SStot;
end
clear SSres SStot

%   Issue warning if R2 is too low
if prctile(R2,90) < 0.9
    tbolus = [];
    return;
end

%   Find best R2 then fit local slope and intercept (this seems to work better than using diff)
ind = find(R2 == max(R2),1);
[~,sl,in] = locallin(tstamp,k(ind,:),9);

%   Find peaks then look for first peak with large amplitude
[pkht,pkloc] = findpeaks(sl,'NPeaks',1,'MinPeakHeight',0.5*max(sl),'MinPeakWidth',3);

%   Use this peak location and slope to regress back to baseline
kbase = median(k(:,1:(pkloc-round(5/dt))),2);
tbolus = (kbase(ind) - in(pkloc))./sl(pkloc);


%   Plotting stuff
if plot_fit
    figure(1021);
    subplot(3,1,1);
    hold off;
    imagesc(log(100*abs(k)./max(abs(k(:)))+1));
    %imagesc(k);
    colormap jet;
    hold on;
    tind = find(abs(tbolus-tstamp) == min(abs(tbolus-tstamp)),1);
    line([tind tind],[1 size(k,1)],'LineWidth',2,'LineStyle','--','Color','k');
    line([1 length(tstamp)],[ind(1) ind(1)],'LineWidth',2,'LineStyle','--','Color','r');
    set(gca,'XTick',[],'YTick',[]);
    ylabel('Readout point');xlabel('Time');
    
    subplot(3,1,3);
    hold off;
    tm = tstamp(tstamp>=tbolus & tstamp<=tstamp(pkloc)+2);
    kslope = polyval([sl(pkloc) in(pkloc)],tm);
    plt2 = plot(tstamp0,k0(ind(1),:),tstamp,kf(ind(1),:),tm,kslope);
    xlabel('Time (s)');ylabel('Signal intensity (arb units)');
    grid on;
    plt2(1).LineWidth = 0.2;
    plt2(3).LineWidth = 4;
    legend('Raw k0 data','Smoothed k0 data','Slope fit');
    hold on;
    line([tbolus tbolus],get(gca,'YLim'),'LineWidth',2,'LineStyle','--','Color','k');
    set(gca,'XLim',[tstamp0(1) tstamp0(end)]);
    
    subplot(3,1,2);
    plt1 = plot(tstamp,sl,tstamp(pkloc),pkht,'o');
    xlabel('Time (s)');ylabel('Local slope (arb units)');
    legend('Local slope','Peak slope at estimated first pass');
    grid on
    set(gca,'XLim',[tstamp0(1) tstamp0(end)]);
    drawnow;
end

%   Plot
if export
    %   Export
    mkdir bolus_data;cd bolus_data;
    dlmwrite('Raw_k0.txt',[(tstamp0-tbolus)' k0(ind(1),:)'],'delimiter',' ','newline','unix');
    dlmwrite('Filtered_k.txt',[(tstamp-tbolus) k(ind(1),:)'],'delimiter',' ','newline','unix');
    dlmwrite('Fit.txt',[tm(:)-tbolus kslope(:)],'delimiter',' ','newline','unix');
    cd ..;
end

end
