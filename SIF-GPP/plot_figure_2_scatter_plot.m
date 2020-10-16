%%% calculate R2 between hotspot and nadir SIF and GPP
%%% written by dalei hao
%%%
clc
clear all;

modelfun = @(b,x) (b(1)*x(:, 1))./(x(:, 1) + b(2));

%% wheat
R2s = zeros(1,4);
P_values = zeros(1, 4);

%% load data
load('../results/sif_gpp_wheat.mat');
FarSIFs_1 = FarSIFs(:,19);
nadirFarSIFs_rs_1 = nadirFarSIF_rs(:,19);
nadirFarSIFs_ks_1 = nadirFarSIF_ks(:,19);
totalFarSIFs_fpars_1 = totalFarSIF_fpars(:,19)*pi;
totalFarSIFs_i0s_1 = totalFarSIF_i0s(:,19)*pi;
RedSIFs_1 = RedSIFs(:,19);
nadirRedSIFs_rs_1 = nadirRedSIF_rs(:,19);
nadirRedSIFs_ks_1 = nadirRedSIF_ks(:,19);
totalRedSIFs_fpars_1 = totalRedSIF_fpars(:,19)*pi;
totalRedSIFs_i0s_1 = totalRedSIF_i0s(:,19)*pi;
GPPs_1 = GPPs(:,19);

filters = GPPs_1>0 & FarSIFs_1>0 & nadirFarSIFs_rs_1>0 & nadirFarSIFs_ks_1>0 &  totalFarSIFs_fpars_1>0 & totalFarSIFs_i0s_1>0 & ...
    GPPs_1<60 & FarSIFs_1<4 & nadirFarSIFs_rs_1<4 & nadirFarSIFs_ks_1<4 &  totalFarSIFs_fpars_1<25 & totalFarSIFs_i0s_1<25 & ...
    RedSIFs_1>0 & nadirRedSIFs_rs_1>0 & nadirRedSIFs_rs_1>0 &  totalRedSIFs_fpars_1>0 & totalRedSIFs_i0s_1>0 & ...
    RedSIFs_1<1 & nadirRedSIFs_rs_1<1 & nadirRedSIFs_rs_1<1 &  totalRedSIFs_fpars_1<12 & totalRedSIFs_i0s_1<12;


FarSIFs_1 = FarSIFs_1(filters);
nadirFarSIFs_rs_1 = nadirFarSIFs_rs_1(filters);
nadirFarSIFs_ks_1 = nadirFarSIFs_ks_1(filters);
totalFarSIFs_fpars_1 = totalFarSIFs_fpars_1(filters);
totalFarSIFs_i0s_1 = totalFarSIFs_i0s_1(filters);
RedSIFs_1 = RedSIFs_1(filters);
nadirRedSIFs_rs_1 = nadirRedSIFs_rs_1(filters);
nadirRedSIFs_ks_1 = nadirRedSIFs_ks_1(filters);
totalRedSIFs_fpars_1 = totalRedSIFs_fpars_1(filters);
totalRedSIFs_i0s_1 = totalRedSIFs_i0s_1(filters);
GPPs_1 = GPPs_1(filters);

load('../results/sif_gpp_corn.mat');
FarSIFs_2 = FarSIFs(:,19);
nadirFarSIFs_rs_2 = nadirFarSIF_rs(:,19);
nadirFarSIFs_ks_2 = nadirFarSIF_ks(:,19);
totalFarSIFs_fpars_2 = totalFarSIF_fpars(:,19)*pi;
totalFarSIFs_i0s_2 = totalFarSIF_i0s(:,19)*pi;
RedSIFs_2 = RedSIFs(:,19);
nadirRedSIFs_rs_2 = nadirRedSIF_rs(:,19);
nadirRedSIFs_ks_2 = nadirRedSIF_ks(:,19);
totalRedSIFs_fpars_2 = totalRedSIF_fpars(:,19)*pi;
totalRedSIFs_i0s_2 = totalRedSIF_i0s(:,19)*pi;
GPPs_2 = GPPs(:,19);

filters = GPPs_2>0 & FarSIFs_2>0 & nadirFarSIFs_rs_2>0 & nadirFarSIFs_ks_2>0 &  totalFarSIFs_fpars_2>0 & totalFarSIFs_i0s_2>0 & ...
    GPPs_2<60 & FarSIFs_2<4 & nadirFarSIFs_rs_2<4 & nadirFarSIFs_ks_2<4 &  totalFarSIFs_fpars_2<25 & totalFarSIFs_i0s_2<25 & ...
    RedSIFs_2>0 & nadirRedSIFs_rs_2>0 & nadirRedSIFs_rs_2>0 &  totalRedSIFs_fpars_2>0 & totalRedSIFs_i0s_2>0 & ...
    RedSIFs_2<1 & nadirRedSIFs_rs_2<1 & nadirRedSIFs_rs_2<1 &  totalRedSIFs_fpars_2<12 & totalRedSIFs_i0s_2<12;


FarSIFs_2 = FarSIFs_2(filters);
nadirFarSIFs_rs_2 = nadirFarSIFs_rs_2(filters);
nadirFarSIFs_ks_2 = nadirFarSIFs_ks_2(filters);
totalFarSIFs_fpars_2 = totalFarSIFs_fpars_2(filters);
totalFarSIFs_i0s_2 = totalFarSIFs_i0s_2(filters);
RedSIFs_2 = RedSIFs_2(filters);
nadirRedSIFs_rs_2 = nadirRedSIFs_rs_2(filters);
nadirRedSIFs_ks_2 = nadirRedSIFs_ks_2(filters);
totalRedSIFs_fpars_2 = totalRedSIFs_fpars_2(filters);
totalRedSIFs_i0s_2 = totalRedSIFs_i0s_2(filters);
GPPs_2 = GPPs_2(filters);

%% calculate R2
R2s(1,1) = calculateR2(FarSIFs_1, GPPs_1);

[b1,bint,r,rint,stats] = regress(GPPs_1, [ ones(size(FarSIFs_1,1),1) FarSIFs_1 ]);
P_values(1, 1) = stats(3);

R2s(1,2) = calculateR2(RedSIFs_1, GPPs_1);
[b2,bint,r,rint,stats] = regress(GPPs_1, [ ones(size(RedSIFs_1,1),1) RedSIFs_1 ]);
P_values(1, 2) = stats(3);

R2s(1,3) = calculateR2(FarSIFs_2, GPPs_2);
[b3,bint,r,rint,stats] = regress(GPPs_2, [ ones(size(FarSIFs_2,1),1) FarSIFs_2 ]);
P_values(1, 3) = stats(3);

R2s(1,4) = calculateR2(RedSIFs_2, GPPs_2);
[b4,bint,r,rint,stats] = regress(GPPs_2, [ ones(size(RedSIFs_2,1),1) RedSIFs_2 ]);
P_values(1, 4) = stats(3);

P_values

%% hyper fitting
R2_hys = zeros(1,4);
P_value_hys = zeros(1, 4);


startPoints = [0 0];
mdl1 = HyperbolicFit(FarSIFs_1, GPPs_1,startPoints);
mdl2 = HyperbolicFit(RedSIFs_1, GPPs_1,startPoints);
mdl3 = HyperbolicFit(FarSIFs_2, GPPs_2,startPoints);
mdl4 = HyperbolicFit(RedSIFs_2, GPPs_2,startPoints);

mdl1
mdl2
mdl3
mdl4


R2_hys(1, 1) = mdl1.Rsquared.Ordinary;
R2_hys(1, 2) = mdl2.Rsquared.Ordinary;
R2_hys(1, 3) = mdl3.Rsquared.Ordinary;
R2_hys(1, 4) = mdl4.Rsquared.Ordinary;

% P_value_hys(1, 1) = stats1.rsquare;
% P_value_hys(1, 2) = stats2.rsquare;
% P_value_hys(1, 3) = stats3.rsquare;
% P_value_hys(1, 4) = stats4.rsquare;

%% plot raw-data
figure;
set(gcf,'unit','normalized','position',[0.2,0.0,0.5,1]);
set(gca, 'Position', [0 0 1 1])
%set label
subplot('position',[0.1 0.8 0.18 0.15])
hold on

f1 = scatter( FarSIFs_1, GPPs_1,10,'k','filled');
f1.MarkerFaceAlpha = 1;
text(3*0.6,0.2*60,['R^2=' num2str(R2s(1,1),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'r')
text(3*0.6,0.07*60,['R^2=' num2str(R2_hys(1,1),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'b')

min_value = 0;
max_value = 10; 

plot([min_value max_value], [ b1(1,1) + min_value*b1(2,1) b1(1,1) + max_value*b1(2,1) ], 'r', 'linewidth', 1.5)

coefficients = mdl1.Coefficients{:, 'Estimate'};
xFitted = 0:0.1:10; 
yFitted = modelfun(coefficients, xFitted(:));
plot(xFitted, yFitted, 'b', 'linewidth', 1.5)


axis([0 3 0 60])
box on
ylabel('GPP (umol m^{-2} s^{-1})')
set(gca,'fontsize',8,'fontname','time new roman')
title('Wheat - Far-red SIF')

%%%%%%%%%% add y label
h = text(-0.43*3,0.2*60,'Raw-Data','fontsize',10,'fontname','time new roman','fontweight','bold');
set(h,'Rotation',90);
text(-0.55*3,1.1*60,'a','fontsize',12,'fontname','time new roman','fontweight','bold')

subplot('position',[0.1+0.23 0.8 0.18 0.15])
hold on

f1 = scatter(RedSIFs_1, GPPs_1,10,'k','filled');
f1.MarkerFaceAlpha = 1;
text(1*0.6,0.2*60,['R^2=' num2str(R2s(1,2),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'r')
text(1*0.6,0.07*60,['R^2=' num2str(R2_hys(1,2),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'b')

%text(0.005,0.8*60,'p < 0.001','fontsize',8,'fontname','time new roman')

min_value = 0;
max_value = 10; 
plot([min_value max_value], [ b2(1,1) + min_value*b2(2,1) b2(1,1) + max_value*b2(2,1) ], 'r', 'linewidth', 1.5)

coefficients = mdl2.Coefficients{:, 'Estimate'};
xFitted = 0:0.1:10; 
yFitted = modelfun(coefficients, xFitted(:));
plot(xFitted, yFitted, 'b', 'linewidth', 1.5)

box on
axis([0 1 0 60])
set(gca,'fontsize',8,'fontname','time new roman')
title('Wheat - Red SIF')
% plot3
subplot('position',[0.1+0.23*2 0.8 0.18 0.15])
hold on

f1 = scatter(FarSIFs_2, GPPs_2,10,'k','filled');
f1.MarkerFaceAlpha = 1;

text(3*0.6,0.2*60,['R^2=' num2str(R2s(1,3),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'r')
text(3*0.6,0.07*60,['R^2=' num2str(R2_hys(1,3),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'b')
min_value =0;
max_value = 10; 
plot([min_value max_value], [ b3(1,1) + min_value*b3(2,1) b3(1,1) + max_value*b3(2,1) ], 'r', 'linewidth', 1.5)

coefficients = mdl3.Coefficients{:, 'Estimate'};
xFitted = 0:0.1:10; 
yFitted = modelfun(coefficients, xFitted(:));
plot(xFitted, yFitted, 'b', 'linewidth', 1.5)

axis([0 3 0 60])
set(gca,'fontsize',8,'fontname','time new roman')
box on
title('Corn - Far-red SIF')

% 4
subplot('position',[0.1+0.23*3 0.8 0.18 0.15])
hold on

f1 = scatter(RedSIFs_2, GPPs_2,10,'k','filled');
f1.MarkerFaceAlpha = 1;
text(1*0.6,0.2*60,['R^2=' num2str(R2s(1,4),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'r')
text(1*0.6,0.07*60,['R^2=' num2str(R2_hys(1,4),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'b')
min_value = 0;
max_value = 10; 
plot([min_value max_value], [ b4(1,1) + min_value*b4(2,1) b4(1,1) + max_value*b4(2,1) ], 'r', 'linewidth', 1.5)


coefficients = mdl4.Coefficients{:, 'Estimate'};
xFitted = 0:0.1:10; 
yFitted = modelfun(coefficients, xFitted(:));
plot(xFitted, yFitted, 'b', 'linewidth', 1.5)

axis([0 1 0 60])
box on
set(gca,'fontsize',8,'fontname','time new roman')
title('Corn - Red SIF')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot nadir 1
R2s = zeros(1,5);
P_values = zeros(1, 4);

R2s(1,1) = calculateR2(nadirFarSIFs_rs_1, GPPs_1);

[b1,bint,r,rint,stats] = regress(GPPs_1, [ ones(size(nadirFarSIFs_rs_1,1),1) nadirFarSIFs_rs_1 ]);
P_values(1, 1) = stats(3);

R2s(1,2) = calculateR2(nadirRedSIFs_rs_1, GPPs_1);
[b2,bint,r,rint,stats] = regress(GPPs_1, [ ones(size(nadirRedSIFs_rs_1,1),1) nadirRedSIFs_rs_1 ]);
P_values(1, 2) = stats(3);

R2s(1,3) = calculateR2(nadirFarSIFs_rs_2, GPPs_2);
[b3,bint,r,rint,stats] = regress(GPPs_2, [ ones(size(nadirFarSIFs_rs_2,1),1) nadirFarSIFs_rs_2 ]);
P_values(1, 3) = stats(3);

R2s(1,4) = calculateR2(nadirRedSIFs_rs_2, GPPs_2);
[b4,bint,r,rint,stats] = regress(GPPs_2, [ ones(size(nadirRedSIFs_rs_2,1),1) nadirRedSIFs_rs_2 ]);
P_values(1, 4) = stats(3);

P_values

%% hyper fitting
R2_hys = zeros(1,4);
P_value_hys = zeros(1, 4);


startPoints = [0 0];
mdl1 = HyperbolicFit(nadirFarSIFs_rs_1, GPPs_1,startPoints);
mdl2 = HyperbolicFit(nadirRedSIFs_rs_1, GPPs_1,startPoints);
mdl3 = HyperbolicFit(nadirFarSIFs_rs_2, GPPs_2,startPoints);
mdl4 = HyperbolicFit(nadirRedSIFs_rs_2, GPPs_2,startPoints);

mdl1
mdl2
mdl3
mdl4


R2_hys(1, 1) = mdl1.Rsquared.Ordinary;
R2_hys(1, 2) = mdl2.Rsquared.Ordinary;
R2_hys(1, 3) = mdl3.Rsquared.Ordinary;
R2_hys(1, 4) = mdl4.Rsquared.Ordinary;

%% start plot
subplot('position',[0.1 0.8-0.18 0.18 0.15])
hold on

f1 = scatter( nadirFarSIFs_rs_1, GPPs_1,10,'k','filled');
f1.MarkerFaceAlpha = 1;
text(3*0.6,0.2*60,['R^2=' num2str(R2s(1,1),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'r')
text(3*0.6,0.07*60,['R^2=' num2str(R2_hys(1,1),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'b')
%text(0.05,0.8*60,'p < 0.001','fontsize',8,'fontname','time new roman')
min_value = 0;
max_value = 10; 

plot([min_value max_value], [ b1(1,1) + min_value*b1(2,1) b1(1,1) + max_value*b1(2,1) ], 'r', 'linewidth', 1.5)
coefficients = mdl1.Coefficients{:, 'Estimate'};
xFitted = 0:0.1:10; 
yFitted = modelfun(coefficients, xFitted(:));
plot(xFitted, yFitted, 'b', 'linewidth', 1.5)
ylabel('GPP (umol m^{-2} s^{-1})')

axis([0 3 0 60])
box on
set(gca,'fontsize',8,'fontname','time new roman')

%%%%%%%%%% add y label
h = text(-0.43*3,0.0*60,'Nadir-NIRv/Redv','fontsize',10,'fontname','time new roman','fontweight','bold');
set(h,'Rotation',90);
text(-0.55*3,1.1*60,'b','fontsize',12,'fontname','time new roman','fontweight','bold')

subplot('position',[0.1+0.23 0.8-0.18 0.18 0.15])
hold on

f1 = scatter(nadirRedSIFs_rs_1, GPPs_1,10,'k','filled');
f1.MarkerFaceAlpha = 1;
text(1*0.6,0.2*60,['R^2=' num2str(R2s(1,2),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'r')
text(1*0.6,0.07*60,['R^2=' num2str(R2_hys(1,2),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'b')
%text(0.005,0.8*60,'p < 0.001','fontsize',8,'fontname','time new roman')

min_value = 0;
max_value = 10; 
plot([min_value max_value], [ b2(1,1) + min_value*b2(2,1) b2(1,1) + max_value*b2(2,1) ], 'r', 'linewidth', 1.5)
coefficients = mdl2.Coefficients{:, 'Estimate'};
xFitted = 0:0.1:10; 
yFitted = modelfun(coefficients, xFitted(:));
plot(xFitted, yFitted, 'b', 'linewidth', 1.5)


box on
axis([0 1 0 60])
set(gca,'fontsize',8,'fontname','time new roman')

% plot3
subplot('position',[0.1+0.23*2 0.8-0.18 0.18 0.15])
hold on

f1 = scatter(nadirFarSIFs_rs_2, GPPs_2,10,'k','filled');
f1.MarkerFaceAlpha = 1;

text(3*0.6,0.2*60,['R^2=' num2str(R2s(1,3),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'r')
text(3*0.6,0.07*60,['R^2=' num2str(R2_hys(1,3),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'b')
%text(0.05*3,0.8*60,'p < 0.001','fontsize',8,'fontname','time new roman')
min_value =0;
max_value = 10; 
plot([min_value max_value], [ b3(1,1) + min_value*b3(2,1) b3(1,1) + max_value*b3(2,1) ], 'r', 'linewidth', 1.5)
coefficients = mdl3.Coefficients{:, 'Estimate'};
xFitted = 0:0.1:10; 
yFitted = modelfun(coefficients, xFitted(:));
plot(xFitted, yFitted, 'b', 'linewidth', 1.5)


axis([0 3 0 60])
set(gca,'fontsize',8,'fontname','time new roman')
box on

% 4
subplot('position',[0.1+0.23*3 0.8-0.18 0.18 0.15])
hold on

f1 = scatter(nadirRedSIFs_rs_2, GPPs_2,10,'k','filled');
f1.MarkerFaceAlpha = 1;
text(1*0.6,0.2*60,['R^2=' num2str(R2s(1,4),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'r')
text(1*0.6,0.07*60,['R^2=' num2str(R2_hys(1,4),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'b')
%text(0.05,0.8*60,'p < 0.001','fontsize',8,'fontname','time new roman')
min_value = 0;
max_value = 10; 
plot([min_value max_value], [ b4(1,1) + min_value*b4(2,1) b4(1,1) + max_value*b4(2,1) ], 'r', 'linewidth', 1.5)
coefficients = mdl4.Coefficients{:, 'Estimate'};
xFitted = 0:0.1:10; 
yFitted = modelfun(coefficients, xFitted(:));
plot(xFitted, yFitted, 'b', 'linewidth', 1.5)

axis([0 1 0 60])
box on
set(gca,'fontsize',8,'fontname','time new roman')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot nadir 2
R2s = zeros(1,5);
P_values = zeros(1, 4);

R2s(1,1) = calculateR2(nadirFarSIFs_ks_1, GPPs_1);

[b1,bint,r,rint,stats] = regress(GPPs_1, [ ones(size(nadirFarSIFs_ks_1,1),1) nadirFarSIFs_ks_1 ]);
P_values(1, 1) = stats(3);

R2s(1,2) = calculateR2(nadirRedSIFs_ks_1, GPPs_1);
[b2,bint,r,rint,stats] = regress(GPPs_1, [ ones(size(nadirRedSIFs_ks_1,1),1) nadirRedSIFs_ks_1 ]);
P_values(1, 2) = stats(3);

R2s(1,3) = calculateR2(nadirFarSIFs_ks_2, GPPs_2);
[b3,bint,r,rint,stats] = regress(GPPs_2, [ ones(size(nadirFarSIFs_ks_2,1),1) nadirFarSIFs_ks_2 ]);
P_values(1, 3) = stats(3);

R2s(1,4) = calculateR2(nadirRedSIFs_ks_2, GPPs_2);
[b4,bint,r,rint,stats] = regress(GPPs_2, [ ones(size(nadirRedSIFs_ks_2,1),1) nadirRedSIFs_ks_2 ]);
P_values(1, 4) = stats(3);


P_values

%% hyper fitting
R2_hys = zeros(1,4);
P_value_hys = zeros(1, 4);


startPoints = [0 0];
mdl1 = HyperbolicFit(nadirFarSIFs_ks_1, GPPs_1,startPoints);
mdl2 = HyperbolicFit(nadirRedSIFs_ks_1, GPPs_1,startPoints);
mdl3 = HyperbolicFit(nadirFarSIFs_ks_2, GPPs_2,startPoints);
mdl4 = HyperbolicFit(nadirRedSIFs_ks_2, GPPs_2,startPoints);

mdl1
mdl2
mdl3
mdl4


R2_hys(1, 1) = mdl1.Rsquared.Ordinary;
R2_hys(1, 2) = mdl2.Rsquared.Ordinary;
R2_hys(1, 3) = mdl3.Rsquared.Ordinary;
R2_hys(1, 4) = mdl4.Rsquared.Ordinary;

%% start plot
subplot('position',[0.1 0.8-0.18*2 0.18 0.15])
hold on

f1 = scatter( nadirFarSIFs_ks_1, GPPs_1,10,'k','filled');
f1.MarkerFaceAlpha = 1;
text(3*0.6,0.2*60,['R^2=' num2str(R2s(1,1),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'r')
text(3*0.6,0.07*60,['R^2=' num2str(R2_hys(1,1),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'b')
%text(0.05,0.8*60,'p < 0.001','fontsize',8,'fontname','time new roman')
min_value = 0;
max_value = 10; 

plot([min_value max_value], [ b1(1,1) + min_value*b1(2,1) b1(1,1) + max_value*b1(2,1) ], 'r', 'linewidth', 1.5)
coefficients = mdl1.Coefficients{:, 'Estimate'};
xFitted = 0:0.1:10; 
yFitted = modelfun(coefficients, xFitted(:));
plot(xFitted, yFitted, 'b', 'linewidth', 1.5)
ylabel('GPP (umol m^{-2} s^{-1})')


axis([0 3 0 60])
box on
set(gca,'fontsize',8,'fontname','time new roman')

%%%%%%%%%% add y label
h = text(-0.43*3,0.2*60,'Nadir-KD','fontsize',10,'fontname','time new roman','fontweight','bold');
set(h,'Rotation',90);
text(-0.55*3,1.1*60,'c','fontsize',12,'fontname','time new roman','fontweight','bold')

subplot('position',[0.1+0.23 0.8-0.18*2 0.18 0.15])
hold on

f1 = scatter(nadirRedSIFs_ks_1, GPPs_1,10,'k','filled');
f1.MarkerFaceAlpha = 1;
text(1*0.6,0.2*60,['R^2=' num2str(R2s(1,2),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'r')
text(1*0.6,0.07*60,['R^2=' num2str(R2_hys(1,2),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'b')
%text(0.005,0.8*60,'p < 0.001','fontsize',8,'fontname','time new roman')

min_value = 0;
max_value = 10; 
plot([min_value max_value], [ b2(1,1) + min_value*b2(2,1) b2(1,1) + max_value*b2(2,1) ], 'r', 'linewidth', 1.5)
coefficients = mdl2.Coefficients{:, 'Estimate'};
xFitted = 0:0.1:10; 
yFitted = modelfun(coefficients, xFitted(:));
plot(xFitted, yFitted, 'b', 'linewidth', 1.5)


box on
axis([0 1 0 60])
set(gca,'fontsize',8,'fontname','time new roman')

% plot3
subplot('position',[0.1+0.23*2 0.8-0.18*2 0.18 0.15])
hold on

f1 = scatter(nadirFarSIFs_ks_2, GPPs_2,10,'k','filled');
f1.MarkerFaceAlpha = 1;

text(3*0.6,0.2*60,['R^2=' num2str(R2s(1,3),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'r')
text(3*0.6,0.07*60,['R^2=' num2str(R2_hys(1,3),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'b')
%text(0.05*3,0.8*60,'p < 0.001','fontsize',8,'fontname','time new roman')
min_value =0;
max_value = 10; 
plot([min_value max_value], [ b3(1,1) + min_value*b3(2,1) b3(1,1) + max_value*b3(2,1) ], 'r', 'linewidth', 1.5)
coefficients = mdl3.Coefficients{:, 'Estimate'};
xFitted = 0:0.1:10; 
yFitted = modelfun(coefficients, xFitted(:));
plot(xFitted, yFitted, 'b', 'linewidth', 1.5)


axis([0 3 0 60])
set(gca,'fontsize',8,'fontname','time new roman')
box on

% 4
subplot('position',[0.1+0.23*3 0.8-0.18*2 0.18 0.15])
hold on

f1 = scatter(nadirRedSIFs_ks_2, GPPs_2,10,'k','filled');
f1.MarkerFaceAlpha = 1;
text(1*0.6,0.2*60,['R^2=' num2str(R2s(1,4),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'r')
text(1*0.6,0.07*60,['R^2=' num2str(R2_hys(1,4),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'b')
%text(0.05,0.8*60,'p < 0.001','fontsize',8,'fontname','time new roman')
min_value = 0;
max_value = 10; 
plot([min_value max_value], [ b4(1,1) + min_value*b4(2,1) b4(1,1) + max_value*b4(2,1) ], 'r', 'linewidth', 1.5)
coefficients = mdl4.Coefficients{:, 'Estimate'};
xFitted = 0:0.1:10; 
yFitted = modelfun(coefficients, xFitted(:));
plot(xFitted, yFitted, 'b', 'linewidth', 1.5)


axis([0 1 0 60])
box on
set(gca,'fontsize',8,'fontname','time new roman')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot total 1
R2s = zeros(1,5);
P_values = zeros(1, 4);

R2s(1,1) = calculateR2(totalFarSIFs_fpars_1, GPPs_1);

[b1,bint,r,rint,stats] = regress(GPPs_1, [ ones(size(totalFarSIFs_fpars_1,1),1) totalFarSIFs_fpars_1 ]);
P_values(1, 1) = stats(3);

R2s(1,2) = calculateR2(totalRedSIFs_fpars_1, GPPs_1);
[b2,bint,r,rint,stats] = regress(GPPs_1, [ ones(size(totalRedSIFs_fpars_1,1),1) totalRedSIFs_fpars_1 ]);
P_values(1, 2) = stats(3);

R2s(1,3) = calculateR2(totalFarSIFs_fpars_2, GPPs_2);
[b3,bint,r,rint,stats] = regress(GPPs_2, [ ones(size(totalFarSIFs_fpars_2,1),1) totalFarSIFs_fpars_2 ]);
P_values(1, 3) = stats(3);

R2s(1,4) = calculateR2(totalRedSIFs_fpars_2, GPPs_2);
[b4,bint,r,rint,stats] = regress(GPPs_2, [ ones(size(totalRedSIFs_fpars_2,1),1) totalRedSIFs_fpars_2 ]);
P_values(1, 4) = stats(3);

P_values

%% hyper fitting
R2_hys = zeros(1,4);
P_value_hys = zeros(1, 4);


startPoints = [0 0];
mdl1 = HyperbolicFit(totalFarSIFs_fpars_1, GPPs_1,startPoints);
mdl2 = HyperbolicFit(totalRedSIFs_fpars_1, GPPs_1,startPoints);
mdl3 = HyperbolicFit(totalFarSIFs_fpars_2, GPPs_2,startPoints);
mdl4 = HyperbolicFit(totalRedSIFs_fpars_2, GPPs_2,startPoints);

mdl1
mdl2
mdl3
mdl4


R2_hys(1, 1) = mdl1.Rsquared.Ordinary;
R2_hys(1, 2) = mdl2.Rsquared.Ordinary;
R2_hys(1, 3) = mdl3.Rsquared.Ordinary;
R2_hys(1, 4) = mdl4.Rsquared.Ordinary;

%% start plot
subplot('position',[0.1 0.8-0.18*3 0.18 0.15])
hold on

f1 = scatter( totalFarSIFs_fpars_1, GPPs_1,10,'k','filled');
f1.MarkerFaceAlpha = 1;
text(25*0.6,0.2*60,['R^2=' num2str(R2s(1,1),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'r')
text(25*0.6,0.07*60,['R^2=' num2str(R2_hys(1,1),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'b')
%text(0.05,0.8*60,'p < 0.001','fontsize',8,'fontname','time new roman')
min_value = 0;
max_value = 25; 

plot([min_value max_value], [ b1(1,1) + min_value*b1(2,1) b1(1,1) + max_value*b1(2,1) ], 'r', 'linewidth', 1.5)
coefficients = mdl1.Coefficients{:, 'Estimate'};
xFitted = 0:0.1:25; 
yFitted = modelfun(coefficients, xFitted(:));
plot(xFitted, yFitted, 'b', 'linewidth', 1.5)

ylabel('GPP (umol m^{-2} s^{-1})')

axis([0 25 0 60])
box on
set(gca,'fontsize',8,'fontname','time new roman')

%%%%%%%%%% add y label
h = text(-0.43*25,0.2*60,'Total-FPAR','fontsize',10,'fontname','time new roman','fontweight','bold');
set(h,'Rotation',90);
text(-0.55*25,1.1*60,'d','fontsize',12,'fontname','time new roman','fontweight','bold')

subplot('position',[0.1+0.23 0.8-0.18*3 0.18 0.15])
hold on

f1 = scatter(totalRedSIFs_fpars_1, GPPs_1,10,'k','filled');
f1.MarkerFaceAlpha = 1;
text(10*0.6,0.2*60,['R^2=' num2str(R2s(1,2),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'r')
text(10*0.6,0.07*60,['R^2=' num2str(R2_hys(1,2),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'b')
%text(0.005,0.8*60,'p < 0.001','fontsize',8,'fontname','time new roman')

min_value = 0;
max_value = 10; 
plot([min_value max_value], [ b2(1,1) + min_value*b2(2,1) b2(1,1) + max_value*b2(2,1) ], 'r', 'linewidth', 1.5)
coefficients = mdl2.Coefficients{:, 'Estimate'};
xFitted = 0:0.1:10; 
yFitted = modelfun(coefficients, xFitted(:));
plot(xFitted, yFitted, 'b', 'linewidth', 1.5)


box on
axis([0 10 0 60])
set(gca,'fontsize',8,'fontname','time new roman')

% plot3
subplot('position',[0.1+0.23*2 0.8-0.18*3 0.18 0.15])
hold on

f1 = scatter(totalFarSIFs_fpars_2, GPPs_2,10,'k','filled');
f1.MarkerFaceAlpha = 1;

text(25*0.6,0.2*60,['R^2=' num2str(R2s(1,3),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'r')
text(25*0.6,0.07*60,['R^2=' num2str(R2_hys(1,3),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'b')
%text(0.05*3,0.8*60,'p < 0.001','fontsize',8,'fontname','time new roman')
min_value =0;
max_value = 25; 
plot([min_value max_value], [ b3(1,1) + min_value*b3(2,1) b3(1,1) + max_value*b3(2,1) ], 'r', 'linewidth', 1.5)
coefficients = mdl3.Coefficients{:, 'Estimate'};
xFitted = 0:0.1:25; 
yFitted = modelfun(coefficients, xFitted(:));
plot(xFitted, yFitted, 'b', 'linewidth', 1.5)


axis([0 25 0 60])
set(gca,'fontsize',8,'fontname','time new roman')
box on

% 4
subplot('position',[0.1+0.23*3 0.8-0.18*3 0.18 0.15])
hold on

f1 = scatter(totalRedSIFs_fpars_2, GPPs_2,10,'k','filled');
f1.MarkerFaceAlpha = 1;
text(10*0.6,0.2*60,['R^2=' num2str(R2s(1,4),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'r')
text(10*0.6,0.07*60,['R^2=' num2str(R2_hys(1,4),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'b')
%text(0.05,0.8*60,'p < 0.001','fontsize',8,'fontname','time new roman')
min_value = 0;
max_value = 10; 
plot([min_value max_value], [ b4(1,1) + min_value*b4(2,1) b4(1,1) + max_value*b4(2,1) ], 'r', 'linewidth', 1.5)
coefficients = mdl4.Coefficients{:, 'Estimate'};
xFitted = 0:0.1:10; 
yFitted = modelfun(coefficients, xFitted(:));
plot(xFitted, yFitted, 'b', 'linewidth', 1.5)


axis([0 10 0 60])
box on
set(gca,'fontsize',8,'fontname','time new roman')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot total 2
R2s = zeros(1,5);
P_values = zeros(1, 4);

R2s(1,1) = calculateR2(totalFarSIFs_i0s_1, GPPs_1);

[b1,bint,r,rint,stats] = regress(GPPs_1, [ ones(size(totalFarSIFs_i0s_1,1),1) totalFarSIFs_i0s_1 ]);
P_values(1, 1) = stats(3);

R2s(1,2) = calculateR2(totalRedSIFs_i0s_1, GPPs_1);
[b2,bint,r,rint,stats] = regress(GPPs_1, [ ones(size(totalRedSIFs_i0s_1,1),1) totalRedSIFs_i0s_1 ]);
P_values(1, 2) = stats(3);

R2s(1,3) = calculateR2(totalFarSIFs_i0s_2, GPPs_2);
[b3,bint,r,rint,stats] = regress(GPPs_2, [ ones(size(totalFarSIFs_i0s_2,1),1) totalFarSIFs_i0s_2 ]);
P_values(1, 3) = stats(3);

R2s(1,4) = calculateR2(totalRedSIFs_i0s_2, GPPs_2);
[b4,bint,r,rint,stats] = regress(GPPs_2, [ ones(size(totalRedSIFs_i0s_2,1),1) totalRedSIFs_i0s_2 ]);
P_values(1, 4) = stats(3);

P_values

%% hyper fitting
R2_hys = zeros(1,4);
P_value_hys = zeros(1, 4);


startPoints = [0 0];
mdl1 = HyperbolicFit(totalFarSIFs_i0s_1, GPPs_1,startPoints);
mdl2 = HyperbolicFit(totalRedSIFs_i0s_1, GPPs_1,startPoints);
mdl3 = HyperbolicFit(totalFarSIFs_i0s_2, GPPs_2,startPoints);
mdl4 = HyperbolicFit(totalRedSIFs_i0s_2, GPPs_2,startPoints);

mdl1
mdl2
mdl3
mdl4


R2_hys(1, 1) = mdl1.Rsquared.Ordinary;
R2_hys(1, 2) = mdl2.Rsquared.Ordinary;
R2_hys(1, 3) = mdl3.Rsquared.Ordinary;
R2_hys(1, 4) = mdl4.Rsquared.Ordinary;

subplot('position',[0.1 0.8-0.18*4 0.18 0.15])
hold on

f1 = scatter( totalFarSIFs_i0s_1, GPPs_1,10,'k','filled');
f1.MarkerFaceAlpha = 1;
text(25*0.6,0.2*60,['R^2=' num2str(R2s(1,1),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'r')
text(25*0.6,0.07*60,['R^2=' num2str(R2_hys(1,1),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'b')
%text(0.05,0.8*60,'p < 0.001','fontsize',8,'fontname','time new roman')
min_value = 0;
max_value = 25; 

plot([min_value max_value], [ b1(1,1) + min_value*b1(2,1) b1(1,1) + max_value*b1(2,1) ], 'r', 'linewidth', 1.5)
coefficients = mdl1.Coefficients{:, 'Estimate'};
xFitted = 0:0.1:25; 
yFitted = modelfun(coefficients, xFitted(:));
plot(xFitted, yFitted, 'b', 'linewidth', 1.5)

ylabel('GPP (umol m^{-2} s^{-1})')

axis([0 25 0 60])
box on
set(gca,'fontsize',8,'fontname','time new roman')

%%%%%%%%%% add y label
h = text(-0.43*25,0.2*60,'Total-i_0','fontsize',10,'fontname','time new roman','fontweight','bold');
set(h,'Rotation',90);
text(-0.55*25,1.1*60,'e','fontsize',12,'fontname','time new roman','fontweight','bold')

subplot('position',[0.1+0.23 0.8-0.18*4 0.18 0.15])
hold on

f1 = scatter(totalRedSIFs_i0s_1, GPPs_1,10,'k','filled');
f1.MarkerFaceAlpha = 1;
text(10*0.6,0.2*60,['R^2=' num2str(R2s(1,2),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'r')
text(10*0.6,0.07*60,['R^2=' num2str(R2_hys(1,2),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'b')
%text(0.005,0.8*60,'p < 0.001','fontsize',8,'fontname','time new roman')

min_value = 0;
max_value = 10; 
plot([min_value max_value], [ b2(1,1) + min_value*b2(2,1) b2(1,1) + max_value*b2(2,1) ], 'r', 'linewidth', 1.5)
coefficients = mdl2.Coefficients{:, 'Estimate'};
xFitted = 0:0.1:10; 
yFitted = modelfun(coefficients, xFitted(:));
plot(xFitted, yFitted, 'b', 'linewidth', 1.5)


box on
axis([0 10 0 60])
set(gca,'fontsize',8,'fontname','time new roman')

% plot3
subplot('position',[0.1+0.23*2 0.8-0.18*4 0.18 0.15])
hold on

f1 = scatter(totalFarSIFs_i0s_2, GPPs_2,10,'k','filled');
f1.MarkerFaceAlpha = 1;

text(25*0.6,0.2*60,['R^2=' num2str(R2s(1,3),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'r')
text(25*0.6,0.07*60,['R^2=' num2str(R2_hys(1,3),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'b')
%text(0.05*3,0.8*60,'p < 0.001','fontsize',8,'fontname','time new roman')
min_value =0;
max_value = 25; 
plot([min_value max_value], [ b3(1,1) + min_value*b3(2,1) b3(1,1) + max_value*b3(2,1) ], 'r', 'linewidth', 1.5)
coefficients = mdl3.Coefficients{:, 'Estimate'};
xFitted = 0:0.1:25; 
yFitted = modelfun(coefficients, xFitted(:));
plot(xFitted, yFitted, 'b', 'linewidth', 1.5)


axis([0 25 0 60])
set(gca,'fontsize',8,'fontname','time new roman')
box on

% 4
subplot('position',[0.1+0.23*3 0.8-0.18*4 0.18 0.15])
hold on

f1 = scatter(totalRedSIFs_i0s_2, GPPs_2,10,'k','filled');
f1.MarkerFaceAlpha = 1;
text(10*0.6,0.2*60,['R^2=' num2str(R2s(1,4),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'r')
text(10*0.6,0.07*60,['R^2=' num2str(R2_hys(1,4),'%4.2f')],'fontsize',8,'fontname','time new roman', 'color', 'b')
%text(0.05,0.8*60,'p < 0.001','fontsize',8,'fontname','time new roman')
min_value = 0;
max_value = 10; 
plot([min_value max_value], [ b4(1,1) + min_value*b4(2,1) b4(1,1) + max_value*b4(2,1) ], 'r', 'linewidth', 1.5)
coefficients = mdl4.Coefficients{:, 'Estimate'};
xFitted = 0:0.1:10; 
yFitted = modelfun(coefficients, xFitted(:));
plot(xFitted, yFitted, 'b', 'linewidth', 1.5)

axis([0 10 0 60])
box on
set(gca,'fontsize',8,'fontname','time new roman')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%



[ax1,h1]=suplabel('Observed/nadir/total SIF (mW m^{-2} nm^{-1} (sr^{-1}))');
ax1.Position(2) = ax1.Position(2) - 0.02;
set(h1,'FontSize',10)

% save figure
print(gcf, '-dtiff', '-r600', 'figure_2_scatter.tif')
close all