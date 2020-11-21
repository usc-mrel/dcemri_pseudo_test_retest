clear all, close all


vp = 0.2;
kt = 0.4 / 60;
kep = 0;

% no negative time
t = -0:5:250;
AIF = AIF_parker(t, 0);

C1 = model_standard(1, length(t), vp, kt, kep, AIF, t, 'sum');
intAIF = integrateAIF(t, AIF, 'sum');
C2 = model_patlak(1, length(t), vp, kt, AIF, intAIF);

figure(1),
plot(t, C1), hold on
plot(t, C2)

% w/ negative time
t = -20:5:250;
AIF = AIF_parker(t, 0);

C3 = model_standard(1, length(t), vp, kt, kep, AIF, t, 'sum');
intAIF = integrateAIF(t, AIF, 'sum');
C4 = model_patlak(1, length(t), vp, kt, AIF, intAIF);

figure(1),
plot(t, C3)
plot(t, C4)

legend({'ETK no neg time', 'Patlak no neg time', 'ETK w/ neg time', 'Patlak w/ neg time'})