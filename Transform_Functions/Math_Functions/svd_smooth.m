function [img2,BF] = svd_smooth(img,n,debug,Vin)

if nargin < 3 || isempty(debug)
    debug = 0;
end

[np, nv, ns, nt] = size(img);

%   Check that n<=nt
if n > nt
    n=nt;
end

img2 = reshape(img,[np*nv*ns nt]);

if nargin < 4 || isempty(Vin)
    [~,S,V] = svd(img2,0);
    BF = V(:,1:n);
    if debug
        figure(1);
        subplot(2,1,1);
        semilogy(diag(S),'x-');
        subplot(2,1,2);
        plot(abs(BF));
        legend(num2str((1:n)'));
        drawnow;
    end
    clear S
else
    BF = Vin(:,1:n);
end

wt = pinv(BF)*img2';

img2 = (BF*wt);
img2 = reshape(img2',size(img));
