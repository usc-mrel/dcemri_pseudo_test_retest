function shiftedAIF = shiftAIF(AIF, time, delay)
%function shiftedAIF = shiftAIF(AIF, time, delay)
%   shifts the time function based on Akima interpolation
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

shiftedAIF = interpolateAIF(AIF, time+delay, time);

end

