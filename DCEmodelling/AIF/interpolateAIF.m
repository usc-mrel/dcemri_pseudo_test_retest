function upsampled_AIF = interpolateAIF(AIF, time, upsampled_time)
%upsampled_AIF = interpolateAIF(time, AIF, upsampled_time)
%   Spline interpolation is common way of shifting AIFs, yet the typical C2
%   nature of B-splines cannot handle the AIF onset so well, hence Akima
%   interpolation
%
%
%   References:
%        1. Nadav G, Liberman G, Artzi M, Kiryati N, Bashat D Ben. Optimization of two-compartment-exchange-model analysis for dynamic contrast-enhanced mri incorporating bolus arrival time. J. Magn. Reson. Imaging 2017;45:237?249 doi: 10.1002/jmri.25362.
%        2. Liberman G, Louzoun Y, Artzi M, Nadav G, Ewing JR, Ben Bashat D. DUSTER: Dynamic contrast enhance up-sampled temporal resolution analysis method. Magn. Reson. Imaging 2016;34:442?450 doi: 10.1016/j.mri.2015.12.014.
%
% Yannick Bliesener 2020

input_shape = size(AIF);

nt = length(AIF);
AIF = reshape(AIF, 1, nt);
time = reshape(time, 1, nt);

output_shape = input_shape;
output_shape(input_shape == nt) = length(upsampled_time);

% determine extrapolation value
AIFtail = AIF(end-10:end);
p = polyfit(1:length(AIFtail),reshape(AIFtail, 1, []),1);
AIFtail = (p(1) * (1:length(AIFtail)) + p(2));
stern_extrapolation_value = mean(AIFtail);
bow_extrapolation_value = 0;

% shift AIF
try 
    upsampled_AIF = interp1(time, AIF, upsampled_time, 'makima', nan);
catch e
    % handle older version of Matlab with a similar C1 estimator
    upsampled_AIF = interp1(time, AIF, upsampled_time, 'pchip', nan);
end

% replace the nan values with respective extrapolation value
split_at = floor(length(upsampled_AIF)/2);
extrapolate_ind = isnan(upsampled_AIF) & [true(1, split_at), false(1, length(upsampled_AIF)-split_at)];
upsampled_AIF(extrapolate_ind) = bow_extrapolation_value;
extrapolate_ind = isnan(upsampled_AIF) & [false(1, split_at), true(1, length(upsampled_AIF)-split_at)];
upsampled_AIF(extrapolate_ind) = stern_extrapolation_value;

upsampled_AIF = reshape(upsampled_AIF, output_shape);
end

