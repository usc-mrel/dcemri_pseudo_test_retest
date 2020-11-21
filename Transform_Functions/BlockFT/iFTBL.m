function img = iFTBL(SP,opt)

sz = size(SP{1}.A);

bl = 16;

img{1} = zeros(sz,'like',SP{1}.A);
img{2} = zeros(sz,'like',SP{1}.A);

[X,Y] = ndgrid((-bl/2+0.5):(bl/2-0.5));
R = sqrt(X.^2+Y.^2)/(bl/2);
R = 1./R;
R = fftshift(fftshift(R,1),2);
clear X Y

for shft = 1:2
    for i = 1:bl:sz(1)
    for j = 1:bl:sz(2)
        ind1 = i:i+bl-1;
        ind2 = j:j+bl-1;
        spSM = SP{shft}.A(ind1,ind2);
        
        %   Scale up DC
        %spSM(1) = spSM(1)*10;
        spSM = spSM .* R;
        
        imSM = ifft2(spSM);
        
        img{shft}(ind1,ind2) = imSM;
    end
    end
end

img{2} = circshift(img{2},[bl/2 bl/2]);
img = (img{1} + img{2}) * 0.5;
