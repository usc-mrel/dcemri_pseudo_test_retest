function B = vesselness(img, sigma)
%function B = vesselness(img, sigma)
%   implements the 'vesselness' metric 
%
%
%   References:
%       [1] Frangi AF, Niessen WJ, Vincken KL, Viergever MA. Multiscale vessel enhancement filtering. 1998;1496:130?137 doi: 10.1007/BFb0056195.
%       [2] Chan SLS, Gal Y. Automatic Detection of Arterial Voxels in Dynamic Contrast-Enhanced MR Images of the Brain. In: International Conference on Digital Image Computing Techniques and Applications (DICTA). ; 2012. pp. 1?7. doi: 10.1109/DICTA.2012.6411710.
%
%   Yannick Bliesener 2020
%   Code largely inspired by MATLAB's fibermetric function
%

switch ndims(img)
    case {0, 1, 2}
        % for these cases just use matlabs build in
        B = fibermetric(img, 6*sigma);
        return
    case 3
        % pass
    otherwise
        error('vesselness: Not implemented for dimenions larger than 3!')
end

% Convert the values to double/single.
classOriginalData = class(img); 
switch (classOriginalData)
case 'uint32'
    img = double(img);
case 'double'
otherwise
    img = single(img);
end

objPolarity = 'bright';

% Constant threshold.
alpha  = 0.5;
beta   = 0.5;
inputc = [];

B = zeros(size(img), 'like', img);

for i = 1:size(sigma,1)
    
    c = inputc;
    
    % compute Hessian of image
    H = hessian3D(img, sigma(i,:));
    H = reshape(H, [], 3, 3);
    H = permute(H, [3 2 1]);
    E = nan(3, size(H,3));
    
    % compute the three eigenvalues at each location
    parfor n = 1:size(H,3)
        E(:, n) = eig(H(:,:,n));
    end
    absE = abs(E);
    
    % eig sorts the eigenvalues in ascending order, yet they need to be
    % sorted in ascending order by absolute value
    [absE, I] = sort(absE, 1, 'ascend');
    I = sub2ind(size(I), I, repmat(1:length(I), [3 1]));
    E = E(I);

    % compute the vesselness metric
    Rb = absE(1,:) ./ sqrt(absE(2,:) .* absE(3,:));
    Ra = absE(2,:) ./ absE(3,:);
    Ssquare = sum(absE.^2, 1);

    if isempty(c)
        maxHessianNorm = max(absE(:));
        c = 0.5*maxHessianNorm;
    end

    V = (1 - exp(-(Ra.^2)/(2*alpha^2))) .* exp(-(Rb.^2)/(2*beta^2)) .* (1 - (exp(-Ssquare/(2*c^2))));

    switch (objPolarity)
    case 'bright'
        V((E(2,:) > 0) | (E(3,:) > 0)) = 0;
    case 'dark'               
        V((E(2,:) < 0) | (E(3,:) < 0)) = 0;
    end
    
    V = reshape(V, size(img));

    % Remove NaN values.
    V(~isfinite(V)) = 0;
    B = max(B, V);
end

% Output should always be single.
B = single(B);

function H = hessian3D(A, sigma)
    
    % first order derivatives
    Gx = gaussianGradient(A, sigma(1), 1);
    Gy = gaussianGradient(A, sigma(2), 2);
    Gz = gaussianGradient(A, sigma(3), 3);
    
    
    % second order derivatives
    Gxx = gaussianGradient(Gx, sigma(1), 1);
    Gyy = gaussianGradient(Gy, sigma(2), 2);
    Gzz = gaussianGradient(Gz, sigma(3), 3);
    
    Gxy = gaussianGradient(Gx, sigma(2), 2);
    Gxz = gaussianGradient(Gx, sigma(3), 3);
    Gyz = gaussianGradient(Gy, sigma(3), 3);
    
    H = cat(ndims(A)+1, Gxx, Gxy, Gxz, Gxy, Gyy, Gyz, Gxz, Gyz, Gzz);
    H = reshape(H, [size(A), 3, 3]);
end

function G = gaussianGradient(A, sigma, axis)
    
    % construct derivative kernel
    x = -ceil(3*sigma):ceil(3*sigma);
    kernel = exp(-0.5*(x./sigma).^2);
    kernel = kernel ./ sum(kernel);
    kernel = -x .* kernel ./ sigma^2;
    
    rvec = ones(1,ndims(A));
    rvec(axis) = length(kernel);
    kernel = reshape(kernel, rvec);
    
    % compute derivative
    G = convn(A, kernel, 'same');
end

end

