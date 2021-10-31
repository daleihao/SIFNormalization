clc;
clear all;
close all;
%% define parameters
site_id = 1;
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
   % allreflectance(allreflectance<=0 | allreflectance>0.2) = nan;
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
cvs = nan(bandnum,4);
cvs_n = nan(bandnum, 4);

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
    cvs(i,1) = std_i./nanmean(Refl_m);
    cvs_n(i,1) = std_n_i./nanmean(Refl_m);
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
    cvs(i,2) = std_i./nanmean(Refl_m);
    cvs_n(i,2) = std_n_i./nanmean(Refl_m);
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
    cvs(i,3) = std_i./nanmean(Refl_m);
    cvs_n(i,3) = std_n_i./nanmean(Refl_m);
%     %% Hapke BRDF
%     Refl_m = allreflectance(i,:);
%     Refl_m = Refl_m';
%     filters = ~isnan(Refl_m);
%     
%     beta0 = [0.2 0 0.1];
%     X = [SZAs VZAs RAAs];
%     y = Refl_m;
%     opts = statset('Display','off','TolFun',1e-15);
%     mdl = fitnlm(X,y,@hapke_kernel,beta0,'Options',opts);
%     
%     %% RTLT kernel
%     Refl_e = mdl.Fitted;
%     
%     allreflectance_est_4(i,:) = Refl_e;
%     
%     coefficients = mdl.Coefficients{:, 'Estimate'};
%     
%     Refl_n  = hapke_kernel(coefficients,[SZA_n VZA_n RAA_n]);
%     
%     allreflectance_n_4(i,:) = Refl_m./Refl_e.*Refl_n;
%     Refl_n = allreflectance_n_4(i,:);
    
%     %% RTLT kernel
%     
%     std_i = std(Refl_m);
%     std_n_i = std(Refl_n);
%     R = corrcoef(Refl_e, Refl_m);
%     RMSE = sqrt(mean((Refl_e-Refl_m).^2));
%     R = R(1,2);
%     R2 = R^2;
%     R2s(i,4) = R2;
%     RMSEs(i,4) = RMSE;
%     stds(i,4) = std_i;
%     stds_n(i,4) = std_n_i;
    
end
%% plot
figure;
set(gcf,'unit','normalized','position',[0.1,0.1,0.6,0.8]);
% 1
subplot('Position', [0.11 0.55 0.37 0.4])
hold on
plot(R2s(:,1),'r','linewidth',1)
plot(R2s(:,2),'G','linewidth',1)
plot(R2s(:,3),'b','linewidth',1)
%plot(R2s(:,4),'k')
legend('RTLT','MWM','MRPV','Location','south')
set(gca,'fontsize',13,'linewidth',1)
%xlabel('Wavelength (nm)')
ylabel('R^2')
box
set(gca,'xtick',[11:100:311], 'xticklabel',[],'fontsize',15,'linewidth',1)
text(10,1*0.1,'(a)','fontsize',20)
ylim([0 1])
% 2

subplot('Position', [0.58 0.55 0.37 0.4])
hold on
plot(RMSEs(:,1),'r','linewidth',1)
plot(RMSEs(:,2),'g','linewidth',1)
plot(RMSEs(:,3),'b','linewidth',1)
legend('RTLT','MWM','MRPV','Location','north')
set(gca,'fontsize',13,'linewidth',1)
ylabel('RMSE')
set(gca,'xtick',[11:100:311], 'xticklabel',[],'fontsize',15,'linewidth',1)
box
text(10,2e-3*0.1,'(b)','fontsize',20)

% 3

subplot('Position', [0.11 0.12 0.37 0.4])
hold on
plot(stds(:,1),'k','linewidth',1)
plot(stds_n(:,1),'r','linewidth',1)
plot(stds_n(:,2),'g','linewidth',1)
plot(stds_n(:,3),'b','linewidth',1)
legend('Observed','RTLT','MWM','MRPV','Location','north')

xlabel('Wavelength (nm)')
ylabel('Std')
box
set(gca,'xtick',[11:100:311], 'xticklabel',[500:500:2000],'fontsize',15,'linewidth',1)
text(10,0.015*0.1,'(c)','fontsize',20)
% 4
subplot('Position', [0.58 0.12 0.37 0.4])
hold on
plot(cvs(:,1),'k','linewidth',1)
plot(cvs_n(:,1),'r','linewidth',1)
plot(cvs_n(:,2),'g','linewidth',1)
plot(cvs_n(:,3),'b','linewidth',1)
legend('Observed','RTLT','MWM','MRPV','Location','north')

xlabel('Wavelength (nm)')
ylabel('CV')
box
set(gca,'xtick',[11:100:311], 'xticklabel',[500:500:2000],'fontsize',15,'linewidth',1)
text(10,0.08*0.1,'(d)','fontsize',20)

print(gcf, '-dtiff', '-r300', 'figure/ce_plot_site1.tif')



%% plot scatter
filters = allreflectance>0 & allreflectance_est_1>0 & allreflectance_est_2>0 & allreflectance_est_3>0;

figure;
cmap = colormap('gray');
colormap(flipud(cmap))
set(gcf,'unit','normalized','position',[0.1,0.1,0.8,0.4]);
% 1
subplot('Position', [0.08 0.2 0.25 0.7])
hold on
plot([0 1],[0 1],'r','linewidth',1)

scatplot(allreflectance(filters), allreflectance_est_1(filters));
set(gca,'fontsize',13,'linewidth',1)
title('(a) RTLT')
xlabel('Observed')
ylabel('Modeled')
xlim([0 0.2])
ylim([0 0.2])
colorbar off
box
R = corrcoef(allreflectance(filters), allreflectance_est_1(filters));
R = R(1,2);
R2 = R^2;
RMSE = sqrt(nanmean((allreflectance(filters)- allreflectance_est_1(filters)).^2));
text(0.02,0.19,['R^2 = ' num2str(R2,'%4.2f')],'fontsize',13)
text(0.02,0.17,['RMSE = ' num2str(RMSE,'%4.3f')],'fontsize',13)

% 2

subplot('Position', [0.37 0.2 0.25 0.7])
hold on
plot([0 1],[0 1],'r','linewidth',1)

scatplot(allreflectance(filters), allreflectance_est_2(filters));
set(gca,'fontsize',13,'linewidth',1,'yticklabel',[])
xlabel('Observed')
%ylabel('Modeled')
colorbar off
box
R = corrcoef(allreflectance(filters), allreflectance_est_2(filters));
R = R(1,2);
R2 = R^2;

RMSE = sqrt(nanmean((allreflectance(filters)- allreflectance_est_2(filters)).^2));
text(0.02,0.19,['R^2 = ' num2str(R2,'%4.2f')],'fontsize',13)
text(0.02,0.17,['RMSE = ' num2str(RMSE,'%4.3f')],'fontsize',13)

xlim([0 0.2])
ylim([0 0.2])
title('(b) MWM')

% 3

subplot('Position', [0.65 0.2 0.25 0.7])
hold on
plot([0 1],[0 1],'r','linewidth',1)

scatplot(allreflectance(filters), allreflectance_est_3(filters));
xlabel('Observed')
title('(c) MRPV')
xlim([0 0.2])
ylim([0 0.2])
box
hcb = colorbar;
hcb.Title.FontWeight = 'Bold';
x1=get(gca,'position');
x=get(hcb,'Position');
% x(3)=0.012;
% x(4)=0.23;
x(1)=0.93;
% x(2)=0.67;
set(hcb,'Position',x')
set(gca,'fontsize',13,'linewidth',1,'yticklabel',[])
R = corrcoef(allreflectance(filters), allreflectance_est_3(filters));
R = R(1,2);
R2 = R^2;

RMSE = sqrt(nanmean((allreflectance(filters)- allreflectance_est_3(filters)).^2));
text(0.02,0.19,['R^2 = ' num2str(R2,'%4.2f')],'fontsize',13)
text(0.02,0.17,['RMSE = ' num2str(RMSE,'%4.3f')],'fontsize',13)

print(gcf, '-dtiff', '-r300', 'figure/ce_scatter_site1.tif')


%% plot image

colors = flipud(brewermap(20, 'Spectral'));

figure;
colormap(colors)
set(gcf,'unit','normalized','position',[0.1,0.1,0.8,0.65]);

% 1
subplot('Position', [0.05 0.56 0.45 0.36])
imagesc(allreflectance', [0 0.2])
set(gca, 'xticklabel',[],'fontsize',15,'linewidth',1)
title('(a) Observed')
% 2

subplot('Position', [0.52 0.56 0.45 0.36])
imagesc(allreflectance_n_1', [0 0.2])
set(gca, 'xticklabel',[],'yticklabel',[],'fontsize',15,'linewidth',1)
colorbar
title('(b) RTLT')
% 3

subplot('Position', [0.05 0.13 0.45 0.36])
imagesc(allreflectance_n_2', [0 0.2])
set(gca,'fontsize',15,'linewidth',1,'xtick',[11:100:311], 'xticklabel',[500:500:2000])
xlabel('Wavelength (nm)')
title('(c) MWM')

% 4

subplot('Position', [0.52 0.13 0.45 0.36])
imagesc(allreflectance_n_3', [0 0.2])
colorbar
set(gca, 'yticklabel',[],'xtick',[11:100:311], 'xticklabel',[500:500:2000],'fontsize',15,'linewidth',1)
title('(d) MRPV')
xlabel('Wavelength (nm)')

print(gcf, '-dtiff', '-r300', 'figure/ce_image_site1.tif')

