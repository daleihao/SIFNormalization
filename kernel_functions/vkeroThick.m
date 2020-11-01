function vol_k= vkeroThick(sZenith,vZenith,rAzimuth)
% Origianl volume kernel Li_SparseR (Roujean et al.,1992; Wanner et al.,1995)
% Angle Input is radian angle
% Code by Dalei Hao,2016/09/21
% Email: dalei.hao.93@gmail.com
%--------------------------------------------------------------------------
%% adjust relative azimuth angle (0-pi)
rAzimuth = abs(rAzimuth);
rAzimuth(rAzimuth>=pi)= 2*pi-rAzimuth(rAzimuth>=pi);

%% calculate kernel value
sz = sZenith;
vz = vZenith;
relaz = rAzimuth;
cosvz = cos(vz);
cossz = cos(sz);
cosxi = cossz.*cosvz + sin(sz).*sin(vz).*cos(relaz);
cosxi(cosxi>=1) = 1;
cosxi(cosxi<=-1) = -1;
xi = acos(cosxi);
K_vol =((pi/2 -xi).*cosxi + sin(xi))./(cossz + cosvz) - pi/4;

%% output
vol_k = K_vol;



