function GRview_samples(tres)

if nargin < 1
    tres = 5;
end

%   Read phase encode table
if exist('GRexttab.txt','file')
    petab = dlmread('GRexttab.txt');
else
    petab = dlmread('GRpetable.txt');
end
npe = size(petab,1);

%   Read MR header
pfiles = dir('P*.7');
pfile = GERecon('Pfile.Load',pfiles(1).name);
hdr = GERecon('Pfile.Header');

%   Get acquisition time stamp
TR = hdr.ImageData.tr;
tstamp = 0:TR:TR*(npe-1);

%   Define viewing bins
tres = tres*1000000;
tview = 0:tres:max(tstamp);


for i = 2:length(tview)
    PE = zeros([pfile.yRes pfile.slices]);
    
    ind = find(tstamp>=tview(i-1) & tstamp<tview(i));
    for j = 1:length(ind);
        y = petab(ind(j),1);
        z = petab(ind(j),2);
        PE(y,z) = 1;
    end
    
    imagesc(PE);
    colormap gray;
    axis square;
    pause(0.2);

end

