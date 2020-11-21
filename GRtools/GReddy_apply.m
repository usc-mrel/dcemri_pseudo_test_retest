function k = GReddy_apply(EC,k,petab)

%   Get sizes
[np, nv, ~, ~, nr] = size(k);

%   Convert full k-space to image space in the readout direction
k = iFastFT(fftshiftF(k,1),1,1);

%   Define spatial position axis (this is the same definition used in GReddy2.m
x = linspace(-0.5,0.5,np)';
X = [ones(np,1) x];

%   Convert phase encode tables from [1 to max] to [-1 to 1]
%   Negate since the pulse sequence uses the negative PE (value of 1 in the table) to measure
PEy = petab(:,1);
PEy = -2*((PEy - 1)./(max(PEy(:))-1) - 0.5);
PEz = petab(:,2);
PEz = -2*((PEz - 1)./(max(PEz(:))-1) - 0.5);


%   Loop through phase encodes and apply phase to subsequent TR period
for isamp = 1:(nv-1)
    
    %   Get donor phase encode values
    pey = PEy(isamp);
    pez = PEz(isamp);    
    
    %   Compute phase shifts occuring at the current sample
    %   Negate since we are removing it
    py = [EC.y.B0; EC.y.Gx];
    pz = [EC.z.B0; EC.z.Gx];
    phi = pey*X*py + pez*X*pz;  %   Evaluate the spatial polynomial and scales by relative pe value
    phi = -phi;
    
    %   Apply phase shifts to k-space data
    phi = exp(sqrt(-1)*phi);
    phi = repmat(phi,[1 1 1 1 nr]);
    k(:,isamp+1,:,:,:) = k(:,isamp+1,:,:,:) .* phi;
    
end

%   Return to full k-space
k = fftshiftF(fFastFT(k,1,1),1);


end

