function M = tile_2dwt(wc)
%   Organize 2-D wavelet decomposition into recursive tiles
%
%   Usage:
%     M = tile_2dwt(wc)
%
%   Input:
%     WC: from fwtN
%
%   Output:
%     M: tiled coefficients
%
%   Author:
%     Travis Smith
%     July, 2011
%

if isstruct(wc)

  tmp = wc;
  clear wc;
  wc{1} = tmp;
  clear tmp;

end


Nlevels = length(wc);
S = zeros(Nlevels,2);
for ii=1:Nlevels
  S(ii,:) = size(wc{ii}.AD);
end

sc = S(:,1).*S(:,2) / (S(end,1)*S(end,2));
sc = sqrt(sc);
sc = ones(size(sc));

for ii=Nlevels-1:-1:1
  S(ii,:) = S(ii+1,:)*2;
end 
  
M = threequad(wc{1},S(1,:)) * sc(1);
for ii=2:Nlevels
  M(1:S(ii-1,1),1:S(ii-1,2)) = threequad(wc{ii},S(ii,:)) * sc(ii);
end
M(1:S(end,1),1:S(end,2)) = wc{end}.AA;


end


function z = threequad(wc,N)

Ny = N(1);
Nx = N(2);
z = zeros(N*2);

[ny,nx] = size(wc.AD);
ix = [1:nx];
iy = [1:ny];

cx = floor(Nx/2) - floor(nx/2);
cy = floor(Ny/2) - floor(ny/2);

% quadrant 1
z(iy+cy,ix+cx+Nx) = wc.AD;

% quadrant 3
z(iy+cy+Ny,ix+cx) = wc.DA;

% quadrant 4
z(iy+cy+Ny,ix+cx+Nx) = wc.DD;

end




