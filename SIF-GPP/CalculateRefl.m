function refl = CalculateRefl(SZA, VZA, RAA, Kparam)

%% deg to rad
SZA = abs(SZA)*pi/180;
VZA = abs(VZA)*pi/180;
RAA = abs(RAA)*pi/180;

count = size(SZA, 1);
Kiso = ones(count, 1);
Kvol = vkeroThick(SZA, VZA, RAA);
Kgeo = LiTransit(SZA, VZA, RAA);

refl = [Kiso Kvol Kgeo]*Kparam;

end