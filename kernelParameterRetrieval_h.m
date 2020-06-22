function [Kparms, R2, RMSE] = kernelParameterRetrieval_h(SZA, VZA, RAA, Refl, band)

%% deg to rad
SZA = abs(SZA)*pi/180;
VZA = abs(VZA)*pi/180;
RAA = abs(RAA)*pi/180;

count = size(SZA, 1);
Kiso = ones(count, 1);
Kvol = vkeroThick_Hotspot_C(SZA, VZA, RAA, band);
Kgeo = LiTransit(SZA, VZA, RAA);
Kparms= [Kiso Kvol Kgeo]\Refl;

Refl_e = [Kiso Kvol Kgeo]*Kparms;
R = corrcoef(Refl_e, Refl);
R = R(1,2);
R2 = R^2;

RMSE = sqrt(mean((Refl_e-Refl).^2));
end