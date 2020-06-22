function vol_k= vkeroThick_Hotspot_C(sZenith,vZenith,rAzimuth,band)
% Hotspot-corrected volume kernel Ross_Thick(Roujean,1992,Wanner,1995; Jiao et al.,2016, 2018)
% Angle Input is radian angle
% hotspot parameters from Jiao et al. 2018 (An algorithm for the retrieval of the clumping index (CI) from the MODIS BRDF product using an adjusted version of the kernel-driven BRDF model)
% band: 1 Red; 2 NIR; 3 Blue;
% Code by Dalei Hao,2016/09/21
% Email: dalei.hao.93@gmail.com
%--------------------------------------------------------------------------

%% set hotspot-corrected papameters based on the specific band
if band == 1 % Red
    c1 =  0.7;
    c2 = 3.1;
elseif band == 2 % NIR
    c1 =  0.6;
    c2 = 3.2;
elseif band == 3 % Blue
    c1 =  0.5;
    c2 = 3.2;
end

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
xi = acos(cosxi);

%% calculate hotspot corrected factor
hotspot_correct_factor = (1+c1*exp(xi*180.0/pi*(-1)/c2));

K_vol =((pi/2 -xi).*cosxi + sin(xi)).*hotspot_correct_factor./(cossz + cosvz) - pi/4;

%% output
vol_k = K_vol;




