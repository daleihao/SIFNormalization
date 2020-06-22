function geo_k = LiTransit(sZenith, vZenith, rAzimuth)
% Origianl geo kernel Li_SparseR (Li et al.,1999)
% Angle Input is radian angle
% Code by Dalei Hao,2016/09/21
% Email: dalei.hao.93@gmail.com
% --------------------------------------------------------------------------
%% adjust relative azimuth angle (0-pi)
rAzimuth = abs(rAzimuth);
rAzimuth(rAzimuth>=pi)= 2*pi-rAzimuth(rAzimuth>=pi);

%% calculate kernel value
brratio = 1;
hbratio = 2;
t1 = brratio*tan(sZenith);
theta_ip = atan(t1);
t2 = brratio*tan(vZenith);
theta_vp = atan(t2);
temp1 = cos(theta_ip);
temp2 = cos(theta_vp);
cosxip = temp1.*temp2+sin(theta_ip).*sin(theta_vp).*cos(rAzimuth);
D1 = tan(theta_ip).*tan(theta_ip)+tan(theta_vp).*tan(theta_vp)-2*tan(theta_ip).*tan(theta_vp).*cos(rAzimuth);
D = sqrt(D1);
cost1 = tan(theta_ip).*tan(theta_vp).*sin(rAzimuth);
cost2 = D1+cost1.*cost1;
temp3 = 1./temp1+1./temp2;
cost = hbratio*sqrt(cost2)./temp3;

cost(cost > 1) = 1;

t = acos(cost);
O = (t-sin(t).*cost).*temp3/pi;
B = temp3-O;

k = zeros(size(sZenith, 1),1);
k(B > 2) = (1+cosxip(B > 2))./(temp2(B > 2).*B(B > 2))-2;
k(B<=2) = -B(B<=2)+(1+cosxip(B<=2))./(2*temp2(B<=2));

%% output
geo_k = k;
end