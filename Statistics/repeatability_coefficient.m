function RC = repeatability_coefficient(data)
%function RC = repeatability_coefficient(data)
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
%       RC = repeatability_coefficient(data);
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

sigma_square = var(data, 0, 2);
RC = 1.96 * sqrt( 2 * mean(sigma_square) );

end

