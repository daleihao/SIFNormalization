function refl = RTLT_model(b, x)
% Origianl RTLT kernel (Hapke et al.,1993; Spurr, 2002)
% Angle Input is radian angle
% Code by Dalei Hao,2021/07/21
% Email: dalei.hao.93@gmail.com
% --------------------------------------------------------------------------
%% input parameters
%% input parameters
k1 = b(1);
k2 = b(2);
k3 = b(3);


SZAs = x(:,1);
VZAs = x(:,2);
RAAs = x(:,3);

RAAs(RAAs>pi) = 2*pi - RAAs(RAAs>pi);

count = size(SZAs, 1);
Kiso = ones(count, 1);
Kvol = vkeroThick(SZAs, VZAs, RAAs);
Kgeo = LiTransit(SZAs, VZAs, RAAs);

refl = Kiso.*k1 + Kvol.*k2 +  Kgeo.*k3;

end
