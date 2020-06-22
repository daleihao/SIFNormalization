function geo_k= gkeroSparseR(sZenith, vZenith, rAzimuth)
% Origianl geo kernel Li_SparseR (Roujean et al.,1992; Wanner et al.,1995; Schaaf  et al.,1994)
% Angle Input is radian angle
% Code by Dalei Hao,2016/09/21
% Email: dalei.hao.93@gmail.com
% --------------------------------------------------------------------------

%% adjust relative azimuth angle (0-pi)
rAzimuth = abs(rAzimuth);
rAzimuth(rAzimuth>=pi)= 2*pi-rAzimuth(rAzimuth>=pi);

%% calculate kernel value
hbratio=2;
hrratio=1;
sZenith=atan(hrratio*tan(sZenith));
vZenith=atan(hrratio*tan(vZenith));
coszeta=cos(sZenith).*cos(vZenith)+sin(sZenith).*sin(vZenith).*cos(rAzimuth);
D2=tan(sZenith).^2+tan(vZenith).^2-2*tan(sZenith).*tan(vZenith).*cos(rAzimuth);
cost=hbratio*sqrt(D2+(tan(sZenith).*tan(vZenith).*sin(rAzimuth)).^2)./(sec(sZenith)+sec(vZenith));
cost(cost(:)>1)=1;
cost(cost(:)<0)=0;
t=acos(cost);
overlap=1/pi*(t-sin(t).*cos(t)).*(sec(sZenith)+sec(vZenith));
LiSparseR=overlap-sec(sZenith)-sec(vZenith)+1/2*(1+coszeta).*sec(vZenith).*sec(sZenith);

%% output
geo_k=LiSparseR;



