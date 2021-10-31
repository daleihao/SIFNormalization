function hapke_k = hapke_kernel(b, x)
% Origianl hapke kernel (Hapke et al.,1993; Spurr, 2002)
% Angle Input is radian angle
% Code by Dalei Hao,2021/07/21
% Email: dalei.hao.93@gmail.com
% --------------------------------------------------------------------------
%% input parameters
omega = b(1);
B0 = b(2);
delta = b(3); 
sZenith = x(:,1);
vZenith = x(:,2);
rAzimuth = x(:,3);

%% adjust relative azimuth angle (0-pi)
rAzimuth = abs(rAzimuth);
rAzimuth(rAzimuth>=pi)= 2*pi-rAzimuth(rAzimuth>=pi);

%% calculate kernel value
ui = cos(sZenith);
ur = cos(vZenith);

cosxi = cos(sZenith).*cos(vZenith) + sin(sZenith).*sin(vZenith).*cos(rAzimuth);
cosxi(cosxi>=1) = 1;
cosxi(cosxi<=-1) = -1;

xi = acos(cosxi);

gamma = sqrt(1-omega);
Rir = omega./ (4*(ui+ur));

P = 1-1.0/2*cosxi;
B = B0.*delta./(delta + tan(1.0/2).*xi);
Ti = (1+2.*ui)./(1+2.*ui.*gamma);
Tr = (1+2.*ur)./(1+2.*ur.*gamma);

%% output
hapke_k = Rir.*((1+B).*P + Ti.*Tr - 1);

