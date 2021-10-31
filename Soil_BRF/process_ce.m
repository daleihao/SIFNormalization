clc;
clear all;
close all;
%% define parameters
site_id = 2;
% sun-sensor geometries
angles_site1 = [68:78
    76.58 75.98 71.36 70.14 69.61 68.90 68.33 63.74 62.86 62.05 61.18
    48.27*ones(1,11)
    346.09*ones(1,11)
    73.73 73.05 67.57 66.05 65.37 64.45 63.70 57.18 55.80 54.51 53.07];

angles_site2 = [79:90
    57.24 57.85 58.26 58.71 59.21 63.47 64.09 64.61 65.19 65.78 66.34 71.08
    44.41*ones(1,12)
    87.70*ones(1,12)
    314.28 313.07 312.26 311.40 310.47 303.33 302.37 301.60 300.74 299.90 299.11 292.86];


%% import hyspectral data
data_ce_observations;

%% select site
if(site_id == 1)% site 1
    allreflectance = allreflectance(:,1:11);
    count = 11;
    angles_site = angles_site1;
    allreflectance(allreflectance<=0 | allreflectance>0.2) = nan;
else%% site 2
    allreflectance = allreflectance(:,12:end);
    count = 12;
    angles_site = angles_site2;
end


for i = 1:count
    allreflectance(:,i) = smoothdata(allreflectance(:,i),'sgolay',51);
end


%% kernel driven fitting
SZAs = angles_site(2,:);
VZAs = angles_site(3,:);
RAAs = abs(angles_site(4,:)-angles_site(5,:));

SZAs = abs(SZAs')*pi/180;
VZAs = abs(VZAs')*pi/180;
RAAs = abs(RAAs')*pi/180;

Kiso = ones(count, 1);
Kvol = vkeroThick(SZAs, VZAs, RAAs);
Kgeo = LiTransit(SZAs, VZAs, RAAs);

SZA_n = ones(count,1)*mean(SZAs);
VZA_n = ones(count,1)*mean(VZAs);
RAA_n = ones(count,1)*mean(RAAs);

Kiso_n = ones(count, 1);
Kvol_n = vkeroThick(SZA_n, VZA_n, RAA_n);
Kgeo_n = LiTransit(SZA_n, VZA_n, RAA_n);

K0_n = ones(count, 1);
K1_n = SZA_n.^2 + VZA_n.^2;
K2_n = (SZA_n.^2).*(VZA_n.^2);
K3_n = cos(RAA_n);

%% modified Walthall
K0 = ones(count, 1);
K1 = SZAs.^2 + VZAs.^2;
K2 = (SZAs.^2).*(VZAs.^2);
K3 = cos(RAAs);

allreflectance_est_1 = zeros(size(allreflectance));
allreflectance_est_2 = zeros(size(allreflectance));
allreflectance_est_3 = zeros(size(allreflectance));
allreflectance_est_4 = zeros(size(allreflectance));

allreflectance_n_1 = zeros(size(allreflectance));
allreflectance_n_2 = zeros(size(allreflectance));
allreflectance_n_3 = zeros(size(allreflectance));
allreflectance_n_4 = zeros(size(allreflectance));

bandnum = size(allreflectance, 1);

R2s = nan(bandnum, 4);
RMSEs = nan(bandnum,4);
stds = nan(bandnum,4);
stds_n = nan(bandnum, 4);
for i = 1:bandnum
    
    Refl_m = allreflectance(i,:);
    Refl_m = Refl_m';
    
    %% RTLT kernel
    Kparms= [Kiso Kvol Kgeo]\Refl_m;
    Refl_e = [Kiso Kvol Kgeo]*Kparms;
    
    allreflectance_est_1(i,:) = Refl_e;
    
    Refl_n = [Kiso_n Kvol_n Kgeo_n]*Kparms;
    
    allreflectance_n_1(i,:) = Refl_m./Refl_e.*Refl_n;
    Refl_n = allreflectance_n_1(i,:);
    
    std_i = std(Refl_m);
    std_n_i = std(Refl_n);
    R = corrcoef(Refl_e, Refl_m);
    RMSE = sqrt(mean((Refl_e-Refl_m).^2));
    R = R(1,2);
    R2 = R^2;
    R2s(i,1) = R2;
    RMSEs(i,1) = RMSE;
    stds(i,1) = std_i;
    stds_n(i,1) = std_n_i;
    
    %% modified Walthall
    Kparms= [K0 K1 K2 K3]\Refl_m;
    Refl_e = [K0 K1 K2 K3]*Kparms;
    
    allreflectance_est_2(i,:) = Refl_e;
    
    Refl_n = [K0_n K1_n K2_n K3_n]*Kparms;
    
    allreflectance_n_2(i,:) = Refl_m./Refl_e.*Refl_n;
    Refl_n = allreflectance_n_2(i,:);
    
    std_i = std(Refl_m);
    std_n_i = std(Refl_n);
    R = corrcoef(Refl_e, Refl_m);
    RMSE = sqrt(mean((Refl_e-Refl_m).^2));
    R = R(1,2);
    R2 = R^2;
    R2s(i,2) = R2;
    RMSEs(i,2) = RMSE;
    stds(i,2) = std_i;
    stds_n(i,2) = std_n_i;
    
    %% modified RPV 
    
    Refl_m = allreflectance(i,:);
    Refl_m = Refl_m';
    
    beta0 = [0.1 1 0.1];
    X = [SZAs VZAs RAAs];
    y = Refl_m;
    opts = statset('Display','off','TolFun',1e-15);
    mdl = fitnlm(X,y,@RPV_kernel,beta0,'Options',opts);
    
    Refl_e = mdl.Fitted;
    
    allreflectance_est_3(i,:) = Refl_e;
    
    coefficients = mdl.Coefficients{:, 'Estimate'};
    
    Refl_n  = RPV_kernel(coefficients,[SZA_n VZA_n RAA_n]);
    
    allreflectance_n_3(i,:) = Refl_m./Refl_e.*Refl_n;
    Refl_n = allreflectance_n_3(i,:);
      
    std_i = std(Refl_m);
    std_n_i = std(Refl_n);
    R = corrcoef(Refl_e, Refl_m);
    RMSE = sqrt(mean((Refl_e-Refl_m).^2));
    R = R(1,2);
    R2 = R^2;
    R2s(i, 3) = R2;
    RMSEs(i,3) = RMSE;
    stds(i,3) = std_i;
    stds_n(i,3) = std_n_i;
    
    %% Hapke BRDF
    Refl_m = allreflectance(i,:);
    Refl_m = Refl_m';
    filters = ~isnan(Refl_m);
    
    beta0 = [0.2 0 0.1];
    X = [SZAs VZAs RAAs];
    y = Refl_m;
    opts = statset('Display','off','TolFun',1e-15);
    mdl = fitnlm(X,y,@hapke_kernel,beta0,'Options',opts);
    
    %% RTLT kernel
    Refl_e = mdl.Fitted;
    
    allreflectance_est_4(i,:) = Refl_e;
    
    coefficients = mdl.Coefficients{:, 'Estimate'};
    
    Refl_n  = hapke_kernel(coefficients,[SZA_n VZA_n RAA_n]);
    
    allreflectance_n_4(i,:) = Refl_m./Refl_e.*Refl_n;
    Refl_n = allreflectance_n_4(i,:);
    
    %% RTLT kernel
    
    std_i = std(Refl_m);
    std_n_i = std(Refl_n);
    R = corrcoef(Refl_e, Refl_m);
    RMSE = sqrt(mean((Refl_e-Refl_m).^2));
    R = R(1,2);
    R2 = R^2;
    R2s(i,4) = R2;
    RMSEs(i,4) = RMSE;
    stds(i,4) = std_i;
    stds_n(i,4) = std_n_i;
    
end
%% plot
figure;
subplot(221)
hold on
plot(R2s(:,1),'r')
plot(R2s(:,2),'G')
plot(R2s(:,3),'b')
legend('RTLT','MWM','MRPV')
title('R2');
subplot(222)
hold on
plot(RMSEs(:,1),'r')
plot(RMSEs(:,2),'g')
plot(RMSEs(:,3),'b')
ylim([0 0.01])

title('RMSE');

figure
subplot(211)
hold on
plot(stds(:,1),'r')
plot(stds_n(:,1),'r')
plot(stds_n(:,2),'g')
plot(stds_n(:,3),'b')

title('std of normalized data');
ylim([0 0.01])
subplot(221)
hold on
plot(stds(:,1),'r')
plot(stds_n(:,1),'r')
plot(stds_n(:,2),'g')
plot(stds_n(:,3),'b')

title('std of normalized data');
ylim([0 0.01])

%% image
figure;
colormap jet
subplot(2,2,1)
imagesc(allreflectance)
title('observed data')
caxis([0 0.2])
subplot(2,2,2)
imagesc(allreflectance_n_1)
caxis([0 0.2])
title('estimated data')
subplot(2,2,3)
imagesc(allreflectance_n_2)
title('normalized data')
caxis([0 0.2])
colorbar

%% plot scatter
figure;
subplot(221)
scatplot(allreflectance_est_1(:), allreflectance(:));
title('RTLT')
subplot(222)
scatplot(allreflectance_est_2(:), allreflectance(:));
title('Modified Walthall')
subplot(223)
scatplot(allreflectance_est_3(:), allreflectance(:));
title('Modified RPV')
subplot(224)
scatplot(allreflectance_est_4(:), allreflectance(:));
title('Hapke')
