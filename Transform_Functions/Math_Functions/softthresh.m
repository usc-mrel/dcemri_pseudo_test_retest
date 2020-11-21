function SP = softthresh(SP,th)

%   Get indices to zero and direction
ind0 = abs(SP) < th;
phase = SP./(abs(SP)+eps);

%   Shrink
SP = SP - th.*phase;

%   Zero out
SP(ind0) = 0;
