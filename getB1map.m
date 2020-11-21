function [b1r, B1par] = getB1map(pname) %#ok<*AGROW>

%   Get current directory
cdir = pwd;

%   Get directory names
if nargin < 1 || isempty(pname)
    for i = 1:2
        pname{i} = uigetdir(pwd,'Select B1 mapping images (alpha first, then 2*alpha)');
    end
end

if length(pname) == 2
    %   Get images
    cd(pname{1});
    [~,a1] = dicomr;
    cd(pname{2});
    [~,a2] = dicomr;
    cd(cdir);
    
    % flip = par2.FlipAngle/2;
    flip = 45;
    
    %   Make relative tip angle map
    b1r = real(180/pi * acos(a2./(2*a1 + eps))) ./ flip;
    
    %   Filter
    filt = gauss2(ones(31,31),0.3);filt = filt./sum(filt(:));
    wt = (a1+a2);
    b1r = convn(b1r.*wt,filt,'same')./(convn(wt,filt,'same')+eps);
    clear filt wt
    
    B1par = [];
    
elseif length(pname) == 1
    cd(pname{1});
    [B1par,b1r] = dicomr;
    cd(cdir);
    [~,~,nsB1] = size(b1r);
    F = repmat(gausswin(11),[1 11]) .* repmat(gausswin(11)',[11 1]);
    F = F / sum(F(:));
    b1r = convn(b1r(:,:,1:nsB1/2).*b1r(:,:,nsB1/2+1:end),F,'same') ./ ...
        (eps + convn(b1r(:,:,nsB1/2+1:end),F,'same'));
    b1r = b1r/10 / B1par.FlipAngle;   %   Converts Block-Siegert output to relative scale factor

end

b1r(b1r==0) = 1; %    Typically occurs due to grad warp
