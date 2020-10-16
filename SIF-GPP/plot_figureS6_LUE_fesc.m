%%% calculate R2 between hotspot and nadir SIF and GPP
%%% written by dalei hao
%%%

clc;
clear all;
%% wheat
R2s = zeros(1,4);
P_values = zeros(1, 4);

%% load data
load('../results/lues_fescs_wheat.mat');
nadirFarFescs = nadirFarFescs(:, 1);
nadirRedFescs = nadirRedFescs(:, 1);
nadirFarFescs_i0 = nadirFarFescs_i0(:, 1);
nadirRedFescs_i0 = nadirRedFescs_i0(:, 1);

LUEs = LUEs(:, 1);
%% calculate R2
R2s(1,1) = calculateR2(LUEs, nadirFarFescs);

[b1,bint,r,rint,stats] = regress(nadirFarFescs, [ ones(size(LUEs,1),1) LUEs ]);
P_values(1, 1) = stats(3);

R2s(1,2) = calculateR2(LUEs, nadirRedFescs);
[b2,bint,r,rint,stats] = regress(nadirRedFescs, [ ones(size(LUEs,1),1) LUEs ]);
P_values(1, 2) = stats(3);

R2s(1,3) = calculateR2(LUEs, nadirFarFescs_i0);

[b3,bint,r,rint,stats] = regress(nadirFarFescs_i0, [ ones(size(LUEs,1),1) LUEs ]);
P_values(1, 3) = stats(3);

R2s(1,4) = calculateR2(LUEs, nadirRedFescs_i0);
[b4,bint,r,rint,stats] = regress(nadirRedFescs_i0, [ ones(size(LUEs,1),1) LUEs ]);
P_values(1, 4) = stats(3);


%% plot
figure;
set(gcf,'unit','normalized','position',[0.2,0.2,0.5,0.4]);
set(gca, 'Position', [0 0 1 1])
%set label
subplot('position',[0.1 0.55 0.18 0.35])
hold on
f1 = scatter( LUEs, nadirFarFescs,10,'k','filled')
f1.MarkerFaceAlpha = 0.5;
text(0.001,0.93*0.5 + 0.3,['R^2=' num2str(R2s(1,1),'%4.2f')],'fontsize',8,'fontname','time new roman')
text(0.001,0.8*0.5 + 0.3,'p<0.001','fontsize',8,'fontname','time new roman')
min_value = 0; 
max_value = 1;

plot([min_value max_value], [ b1(1,1) + min_value*b1(2,1) b1(1,1) + max_value*b1(2,1) ], 'r', 'linewidth', 1.5)
axis([0 0.1 0.3 0.8])
box on
set(gca,'fontsize',8,'fontname','time new roman')
title('Far-red fesc - FPAR')

%%%%%%%%%% add y label
h = text(-0.045,0.5,'Wheat','fontsize',10,'fontname','time new roman','fontweight','bold');
set(h,'Rotation',90);
text(-0.05,0.9,'a','fontsize',12,'fontname','time new roman','fontweight','bold')

% plot 2
subplot('position',[0.1+0.23 0.55 0.18 0.35])
hold on

f1 = scatter(LUEs, nadirRedFescs,10,'k','filled')
f1.MarkerFaceAlpha = 0.5;
text(0.001,0.93*0.3+0.2,['R^2=' num2str(R2s(1,2),'%4.2f')],'fontsize',8,'fontname','time new roman')
text(0.001,0.8*0.3+0.2,'p<0.05','fontsize',8,'fontname','time new roman')

min_value = 0;
max_value =1;
plot([min_value max_value], [ b2(1,1) + min_value*b2(2,1) b2(1,1) + max_value*b2(2,1) ], 'r', 'linewidth', 1.5)
box on
axis([0 0.1 0.2 0.5])
set(gca,'fontsize',8,'fontname','time new roman')
title('Red fesc - FPAR')

%% plot 3
subplot('position',[0.1+0.23*2 0.55 0.18 0.35])
hold on
f1 = scatter( LUEs, nadirFarFescs_i0,10,'k','filled')
f1.MarkerFaceAlpha = 0.5;
text(0.001,0.93*0.5 + 0.3,['R^2=' num2str(R2s(1,3),'%4.2f')],'fontsize',8,'fontname','time new roman')
text(0.001,0.8*0.5 + 0.3,'p<0.001','fontsize',8,'fontname','time new roman')
min_value = 0; 
max_value = 1;

plot([min_value max_value], [ b3(1,1) + min_value*b3(2,1) b3(1,1) + max_value*b3(2,1) ], 'r', 'linewidth', 1.5)
axis([0 0.1 0.3 0.8])
box on
set(gca,'fontsize',8,'fontname','time new roman')
title('Far-red fesc - i_0')

% plot 4
subplot('position',[0.1+0.23*3 0.55 0.18 0.35])
hold on

f1 = scatter(LUEs, nadirRedFescs_i0,10,'k','filled')
f1.MarkerFaceAlpha = 0.5;
text(0.001,0.93*0.3+0.2,['R^2=' num2str(R2s(1,4),'%4.2f')],'fontsize',8,'fontname','time new roman')
text(0.001,0.8*0.3+0.2,'p<0.001','fontsize',8,'fontname','time new roman')

min_value = 0;
max_value =1;
plot([min_value max_value], [ b4(1,1) + min_value*b4(2,1) b4(1,1) + max_value*b4(2,1) ], 'r', 'linewidth', 1.5)
box on
axis([0 0.1 0.2 0.5])
set(gca,'fontsize',8,'fontname','time new roman')
title('Red fesc - i_0')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% corn
R2s = zeros(1,5);
P_values = zeros(1, 4);

% disp(1)
%% load data
load('../results/lues_fescs_corn.mat');
nadirFarFescs = nadirFarFescs(:, 1);
nadirRedFescs = nadirRedFescs(:, 1);
nadirFarFescs_i0 = nadirFarFescs_i0(:, 1);
nadirRedFescs_i0 = nadirRedFescs_i0(:, 1);

LUEs = LUEs(:, 1);
%% calculate R2
R2s(1,1) = calculateR2(LUEs, nadirFarFescs);

[b1,bint,r,rint,stats] = regress(nadirFarFescs, [ ones(size(LUEs,1),1) LUEs ]);
P_values(1, 1) = stats(3);

R2s(1,2) = calculateR2(LUEs, nadirRedFescs);
[b2,bint,r,rint,stats] = regress(nadirRedFescs, [ ones(size(LUEs,1),1) LUEs ]);
P_values(1, 2) = stats(3);

R2s(1,3) = calculateR2(LUEs, nadirFarFescs_i0);

[b3,bint,r,rint,stats] = regress(nadirFarFescs_i0, [ ones(size(LUEs,1),1) LUEs ]);
P_values(1, 3) = stats(3);

R2s(1,4) = calculateR2(LUEs, nadirRedFescs_i0);
[b4,bint,r,rint,stats] = regress(nadirRedFescs_i0, [ ones(size(LUEs,1),1) LUEs ]);
P_values(1, 4) = stats(3);


%% plot
%set label
subplot('position',[0.1 0.12 0.18 0.35])
hold on
f1 = scatter( LUEs, nadirFarFescs,10,'k','filled')
f1.MarkerFaceAlpha = 0.5;
text(0.06,0.2*0.5 + 0.3,['R^2=' num2str(R2s(1,1),'%4.2f')],'fontsize',8,'fontname','time new roman')
text(0.06,0.07*0.5 + 0.3,'p<0.001','fontsize',8,'fontname','time new roman')
min_value = 0; 
max_value = 1;

plot([min_value max_value], [ b1(1,1) + min_value*b1(2,1) b1(1,1) + max_value*b1(2,1) ], 'r', 'linewidth', 1.5)
axis([0 0.1 0.3 0.8])
box on
set(gca,'fontsize',8,'fontname','time new roman')


%%%%%%%%%% add y label
h = text(-0.045,0.5,'Corn','fontsize',10,'fontname','time new roman','fontweight','bold');
set(h,'Rotation',90);
text(-0.05,0.9,'b','fontsize',12,'fontname','time new roman','fontweight','bold')

% plot 2
subplot('position',[0.1+0.23 0.12 0.18 0.35])
hold on

f1 = scatter(LUEs, nadirRedFescs,10,'k','filled')
f1.MarkerFaceAlpha = 0.5;
text(0.06,0.2*0.3 + 0.2,['R^2=' num2str(R2s(1,2),'%4.2f')],'fontsize',8,'fontname','time new roman')
text(0.06,0.07*0.3 + 0.2,'p<0.05','fontsize',8,'fontname','time new roman')

min_value = 0;
max_value =1;
plot([min_value max_value], [ b2(1,1) + min_value*b2(2,1) b2(1,1) + max_value*b2(2,1) ], 'r', 'linewidth', 1.5)
box on
axis([0 0.1 0.2 0.5])
set(gca,'fontsize',8,'fontname','time new roman')


%% plot 3
subplot('position',[0.1+0.23*2 0.12 0.18 0.35])
hold on
f1 = scatter( LUEs, nadirFarFescs_i0,10,'k','filled')
f1.MarkerFaceAlpha = 0.5;
text(0.06,0.2*0.5 + 0.3,['R^2=' num2str(R2s(1,3),'%4.2f')],'fontsize',8,'fontname','time new roman')
text(0.06,0.07*0.5 + 0.3,'p<0.05','fontsize',8,'fontname','time new roman')
min_value = 0; 
max_value = 1;

plot([min_value max_value], [ b3(1,1) + min_value*b3(2,1) b3(1,1) + max_value*b3(2,1) ], 'r', 'linewidth', 1.5)
axis([0 0.1 0.3 0.8])
box on
set(gca,'fontsize',8,'fontname','time new roman')


% plot 4
subplot('position',[0.1+0.23*3 0.12 0.18 0.35])
hold on

f1 = scatter(LUEs, nadirRedFescs_i0,10,'k','filled')
f1.MarkerFaceAlpha = 0.5;
text(0.06,0.2*0.3 + 0.2,['R^2=' num2str(R2s(1,4),'%4.2f')],'fontsize',8,'fontname','time new roman')
text(0.06,0.07*0.3+0.2,'p<0.001','fontsize',8,'fontname','time new roman')

min_value = 0;
max_value =1;
plot([min_value max_value], [ b4(1,1) + min_value*b4(2,1) b4(1,1) + max_value*b4(2,1) ], 'r', 'linewidth', 1.5)
box on
axis([0 0.1 0.2 0.5])
set(gca,'fontsize',8,'fontname','time new roman')


[ax1,h1]=suplabel('LUE');
[ax2,h2]=suplabel('fesc','y');
ax2.Position(1) = ax2.Position(1) - 0.02;
set(h1,'FontSize',10)
set(h2,'FontSize',10)

% save figure
print(gcf, '-dtiff', '-r600', 'figure_S6_LUE_fesc.tif')
close all