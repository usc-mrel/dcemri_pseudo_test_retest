function cov = coefficient_of_variation(data)
%function cov = coefficient_of_variation(data)
%   the coefficient of variation
%
%   Input:
%       data    numeric array 
%                   rows: different subjects/objects
%                   cols: repeat measurements
%
%   Output:
%       the desired statistic
%
%   Example:
%       % 10 Gaussian people with three (perfect) measurements each
%       data = repmat(randn(10, 1), [1 3]);
%       cov = coefficient_of_variation(data);
%
%
% References:
%   1. Obuchowski NA, Reeves AP, Huang EP, et al. Quantitative imaging biomarkers: A review of statistical methods for computer algorithm comparisons. Stat. Methods Med. Res. 2015;24:68?106 doi: 10.1177/0962280214537390.
%       Page 30: Eq. 24 and 25
%   2. Shukla-Dave A, Obuchowski NA, Chenevert TL, et al. Quantitative imaging biomarkers alliance (QIBA) recommendations for improved precision of DWI and DCE-MRI derived biomarkers in multicenter oncology trials. J. Magn. Reson. Imaging 2019;49:e101?e121 doi: 10.1002/jmri.26518.
%
%
% Yannick Bliesener 2020
%

mu = mean(data, 2); % within subject mean
sigma_square = var(data, 0, 2); % within subject variance

% cov = sqrt( mean(sigma_square) ) ./ (mean(mu) + eps); % COV as in Ref 1, Eq. 25 and the following
cov = sqrt( mean( sigma_square ./ (mu.^2 + eps))); % COV as in Ref 2, Table 2

end

