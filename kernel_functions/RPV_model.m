function Refl = RPV_model(parms, angles)
% RPV models fro Rahman et al. (1993)
% Angle Input is radian angle
% Code by Dalei Hao,2021/07/21
% Email: dalei.hao.93@gmail.com
% --------------------------------------------------------------------------
%% input parameters
p0 = parms(1);
k = parms(2);
theta = parms(3);
pc = p0;

sZenith = angles(:,1);
vZenith = angles(:,2);
rAzimuth = angles(:,3);

%% adjust relative azimuth angle (0-pi)
rAzimuth = abs(rAzimuth);
rAzimuth(rAzimuth>=pi)= 2*pi-rAzimuth(rAzimuth>=pi);

%% calculate reflectance value
ui = cos(sZenith);
ur = cos(vZenith);

cosxi = cos(sZenith).*cos(vZenith) + sin(sZenith).*sin(vZenith).*cos(rAzimuth);
cosxi(cosxi>=1) = 1;
cosxi(cosxi<=-1) = -1;
%xi = acos(cosxi);

M1 = power(ui,k-1).*power(ur,k-1)./power(ui+ur,1-k);

FHG = (1-theta.^2)./power(1+theta.^2+2*theta.*cosxi, 1.5);
% FHG = exp(bm*cosxi); //modified RPV

G = power(power(tan(sZenith),2)+power(tan(vZenith),2)-2*tan(sZenith).*tan(vZenith).*cos(rAzimuth),0.5);
H = 1 + (1-pc)./(1+G);
%% output
Refl = p0.*M1.*FHG.*H;

