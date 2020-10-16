function c_factor = CalculateCorrectionFactor(SZA, VZA, RAA, Kparam)

%% deg to rad
SZA = abs(SZA)*pi/180;
VZA = abs(VZA)*pi/180;
RAA = abs(RAA)*pi/180;

count = size(SZA, 1);
Kiso = ones(count, 1);
Kvol = vkeroThick(SZA, VZA, RAA);
Kgeo = LiTransit(SZA, VZA, RAA);
Kvol_0 = vkeroThick(SZA, zeros(count, 1),  zeros(count, 1));
Kgeo_0 = LiTransit(SZA,  zeros(count, 1),  zeros(count, 1));

c_factor = ([Kiso Kvol_0 Kgeo_0]*Kparam) ./ ([Kiso Kvol Kgeo]*Kparam);

end