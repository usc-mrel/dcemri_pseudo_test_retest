function SP = fFTBL(img,opt)

sz = size(img);

bl = 16;

[X,Y] = ndgrid((-bl/2+0.5):(bl/2-0.5));
R = sqrt(X.^2+Y.^2)/(bl/2);
R = fftshift(fftshift(R,1),2);
clear X Y


SP{1}.A = zeros(sz,'like',img);
SP{2}.A = zeros(sz,'like',img);

for shft = 1:2
    if shft == 2
        img = circshift(img,[-bl/2 -bl/2]);
    end
    
    for i = 1:bl:sz(1)
    for j = 1:bl:sz(2)
        ind1 = i:i+bl-1;
        ind2 = j:j+bl-1;
        imSM = img(ind1,ind2);
        
        spSM = fft2(imSM);
        
        
        %   Scale down DC
        %spSM(1) = spSM(1)/10;
        spSM = spSM .* R;
    
        SP{shft}.A(ind1,ind2) = spSM;
        
    end
	end
end
