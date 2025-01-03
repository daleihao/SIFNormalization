function refl = RTLT_model(b, x)
% Origianl RTLT kernel (Hapke et al.,1993; Spurr, 2002)
% Angle Input is radian angle
% Code by Dalei Hao,2021/07/21
% Email: dalei.hao.93@gmail.com
% --------------------------------------------------------------------------
%% input parameters
%% input parameters
k0 = b(1);
k1 = b(2);
k2 = b(3);
k3 = b(4);


SZAs = x(:,1);
VZAs = x(:,2);
RAAs = x(:,3);


RAAs(RAAs>pi) = 2*pi - RAAs(RAAs>pi);

count = size(SZAs, 1);
K0 = ones(count, 1);
K1 = SZAs.^2 + VZAs.^2;
K2 = (SZAs.^2).*(VZAs.^2);
K3 = cos(RAAs);

refl = K0.*k0 + K1.*k1 + K2.*k2 +  K3.*k3;

end
