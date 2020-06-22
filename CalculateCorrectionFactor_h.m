function c_factor = CalculateCorrectionFactor_h(SZA, VZA, RAA, Kparam, band)

%% deg to rad
SZA = abs(SZA)*pi/180;
VZA = abs(VZA)*pi/180;
RAA = abs(RAA)*pi/180;

count = size(SZA, 1);
Kiso = ones(count, 1);
Kvol = vkeroThick(SZA, VZA, RAA);
Kgeo = LiTransit(SZA, VZA, RAA);
Kvol_0 = vkeroThick_Hotspot_C(SZA, zeros(count, 1),  zeros(count, 1), band);
Kgeo_0 = LiTransit(SZA,  zeros(count, 1),  zeros(count, 1));

c_factor = ([Kiso Kvol_0 Kgeo_0]*Kparam) ./ ([Kiso Kvol Kgeo]*Kparam);

end