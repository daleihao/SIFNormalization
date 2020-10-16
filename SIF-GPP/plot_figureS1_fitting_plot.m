%%% calculate R2 between hotspot and nadir SIF and GPP
%%% written by dalei hao
%%%


%% wheat
R2s = zeros(1,4);
P_values = zeros(1, 4);

%% load data
load('../results/fitting_obs_wheat.mat');
fitNIRs = fitNIRs(:);
NIRs = NIRs(:);
fitReds = fitReds(:);
Reds = Reds(:);
fitFarSIFs = fitFarSIFs(:);
FarSIFs = FarSIFs(:);
fitRedSIFs = fitRedSIFs(:);
RedSIFs = RedSIFs(:);
%% calculate R2
R2s(1,1) = calculateR2(fitNIRs, NIRs);

[b1,bint,r,rint,stats] = regress(fitNIRs, [ ones(size(NIRs,1),1) NIRs ]);
P_values(1, 1) = stats(3);

R2s(1,2) = calculateR2(fitReds, Reds);
[b2,bint,r,rint,stats] = regress(fitReds, [ ones(size(Reds,1),1) Reds ]);
P_values(1, 2) = stats(3);

R2s(1,3) = calculateR2(fitFarSIFs, FarSIFs);
[b3,bint,r,rint,stats] = regress(fitFarSIFs, [ ones(size(FarSIFs,1),1) FarSIFs ]);
P_values(1, 3) = stats(3);

R2s(1,4) = calculateR2(fitRedSIFs, RedSIFs);
[b4,bint,r,rint,stats] = regress(fitRedSIFs, [ ones(size(RedSIFs,1),1) RedSIFs ]);
P_values(1, 4) = stats(3);



%% plot
figure;
set(gcf,'unit','normalized','position',[0.2,0.2,0.5,0.4]);
set(gca, 'Position', [0 0 1 1])
%set label
subplot('position',[0.1 0.6 0.18 0.35])
hold on
plot([0 1], [ 0 1 ], 'k', 'linewidth', 0.5);
f1 = scatter( NIRs, fitNIRs,5,'k','filled')
f1.MarkerFaceAlpha = 0.2;
text(0.05,0.93,['R^2=' num2str(R2s(1,1),'%4.2f')],'fontsize',8,'fontname','time new roman')
text(0.05,0.8,'p < 0.001','fontsize',8,'fontname','time new roman')
min_value = 0; 
max_value = max(NIRs)+ 0.1;

plot([min_value max_value], [ b1(1,1) + min_value*b1(2,1) b1(1,1) + max_value*b1(2,1) ], 'r', 'linewidth', 1.5)
axis([0 1 0 1])
box on
set(gca,'fontsize',8,'fontname','time new roman')
title('NIR Reflectance')

%%%%%%%%%% add y label
h = text(-0.45,0.3,'Wheat','fontsize',10,'fontname','time new roman','fontweight','bold');
set(h,'Rotation',90);
text(-0.55,1,'a','fontsize',12,'fontname','time new roman','fontweight','bold')

subplot('position',[0.1+0.23 0.6 0.18 0.35])
hold on
plot([0 1], [ 0 1 ], 'k', 'linewidth', 0.5);
f1 = scatter(Reds, fitReds,5,'k','filled')
f1.MarkerFaceAlpha = 0.2;
text(0.005,0.093,['R^2=' num2str(R2s(1,2),'%4.2f')],'fontsize',8,'fontname','time new roman')
text(0.005,0.08,'p < 0.001','fontsize',8,'fontname','time new roman')

min_value = max([min(Reds)-0.1, 0]);
max_value = max(Reds)+ 0.1;
plot([min_value max_value], [ b2(1,1) + min_value*b2(2,1) b2(1,1) + max_value*b2(2,1) ], 'r', 'linewidth', 1.5)
box on
axis([0 0.1 0 0.1])
set(gca,'fontsize',8,'fontname','time new roman')
title('Red Reflectance')

% plot3
subplot('position',[0.1+0.23*2 0.6 0.18 0.35])
hold on
plot([0 3], [ 0 3 ], 'k', 'linewidth', 0.5);
f1 = scatter(FarSIFs, fitFarSIFs,5,'k','filled')
f1.MarkerFaceAlpha = 0.2;

text(0.05*3,0.93*3,['R^2=' num2str(R2s(1,3),'%4.2f')],'fontsize',8,'fontname','time new roman')
text(0.05*3,0.8*3,'p < 0.001','fontsize',8,'fontname','time new roman')
min_value = max([min(FarSIFs)-0.1, 0]);
max_value = max(FarSIFs)+ 0.1;
plot([min_value max_value], [ b3(1,1) + min_value*b3(2,1) b3(1,1) + max_value*b3(2,1) ], 'r', 'linewidth', 1.5)
axis([0 3 0 3])
set(gca,'fontsize',8,'fontname','time new roman')
box on
title('Far-red SIF')

% 4
subplot('position',[0.1+0.23*3 0.6 0.18 0.35])
hold on
plot([0 1], [ 0 1 ], 'k', 'linewidth', 0.5);
f1 = scatter(RedSIFs, fitRedSIFs,5,'k','filled')
f1.MarkerFaceAlpha = 0.2;
text(0.05,0.93,['R^2=' num2str(R2s(1,4),'%4.2f')],'fontsize',8,'fontname','time new roman')
text(0.05,0.8,'p < 0.001','fontsize',8,'fontname','time new roman')
min_value = max([min(RedSIFs)-0.1, 0]);
max_value = max(RedSIFs)+ 0.1;
plot([min_value max_value], [ b4(1,1) + min_value*b4(2,1) b4(1,1) + max_value*b4(2,1) ], 'r', 'linewidth', 1.5)
axis([0 1 0 1])
box on
set(gca,'fontsize',8,'fontname','time new roman')
title('Red SIF')

%suptitle(['Hourly'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% corn
R2s = zeros(1,5);
P_values = zeros(1, 4);

% disp(1)
%% load data
load('../results/fitting_obs_corn.mat');
fitNIRs = fitNIRs(:);
NIRs = NIRs(:);
fitReds = fitReds(:);
Reds = Reds(:);
fitFarSIFs = fitFarSIFs(:);
FarSIFs = FarSIFs(:);
fitRedSIFs = fitRedSIFs(:);
RedSIFs = RedSIFs(:);
%% calculate R2
R2s(1,1) = calculateR2(fitNIRs, NIRs);

[b1,bint,r,rint,stats] = regress(fitNIRs, [ ones(size(NIRs,1),1) NIRs ]);
P_values(1, 1) = stats(3);

R2s(1,2) = calculateR2(fitReds, Reds);
[b2,bint,r,rint,stats] = regress(fitReds, [ ones(size(Reds,1),1) Reds ]);
P_values(1, 2) = stats(3);

R2s(1,3) = calculateR2(fitFarSIFs, FarSIFs);
[b3,bint,r,rint,stats] = regress(fitFarSIFs, [ ones(size(FarSIFs,1),1) FarSIFs ]);
P_values(1, 3) = stats(3);

R2s(1,4) = calculateR2(fitRedSIFs, RedSIFs);
[b4,bint,r,rint,stats] = regress(fitRedSIFs, [ ones(size(RedSIFs,1),1) RedSIFs ]);
P_values(1, 4) = stats(3);



%% plot
%set label
subplot('position',[0.1 0.15 0.18 0.35])
hold on
plot([0 1], [ 0 1 ], 'k', 'linewidth', 0.5);
f1 = scatter( NIRs, fitNIRs,5,'k','filled')
f1.MarkerFaceAlpha = 0.2;
text(0.05,0.93,['R^2=' num2str(R2s(1,1),'%4.2f')],'fontsize',8,'fontname','time new roman')
text(0.05,0.8,'p < 0.001','fontsize',8,'fontname','time new roman')
min_value = 0; 
max_value = max(NIRs)+ 0.1;

plot([min_value max_value], [ b1(1,1) + min_value*b1(2,1) b1(1,1) + max_value*b1(2,1) ], 'r', 'linewidth', 1.5)
axis([0 1 0 1])
box on
set(gca,'fontsize',8,'fontname','time new roman')
%ylabel('Fitting Values')

%%%%%%%%%% add y label
h = text(-0.45,0.3,'Corn','fontsize',10,'fontname','time new roman','fontweight','bold');
set(h,'Rotation',90);
text(-0.55,1,'b','fontsize',12,'fontname','time new roman','fontweight','bold')

subplot('position',[0.1+0.23 0.15 0.18 0.35])
hold on
plot([0 1], [ 0 1 ], 'k', 'linewidth', 0.5);
f1 = scatter(Reds, fitReds,5,'k','filled')
f1.MarkerFaceAlpha = 0.2;
text(0.005,0.093,['R^2=' num2str(R2s(1,2),'%4.2f')],'fontsize',8,'fontname','time new roman')
text(0.005,0.08,'p < 0.001','fontsize',8,'fontname','time new roman')

min_value = max([min(Reds)-0.1, 0]);
max_value = max(Reds)+ 0.1;
plot([min_value max_value], [ b2(1,1) + min_value*b2(2,1) b2(1,1) + max_value*b2(2,1) ], 'r', 'linewidth', 1.5)
box on
axis([0 0.1 0 0.1])
set(gca,'fontsize',8,'fontname','time new roman')

% plot3
subplot('position',[0.1+0.23*2 0.15 0.18 0.35])
hold on
plot([0 3], [ 0 3 ], 'k', 'linewidth', 0.5);
f1 = scatter(FarSIFs, fitFarSIFs,5,'k','filled')
f1.MarkerFaceAlpha = 0.2;

text(0.05*3,0.93*3,['R^2=' num2str(R2s(1,3),'%4.2f')],'fontsize',8,'fontname','time new roman')
text(0.05*3,0.8*3,'p < 0.001','fontsize',8,'fontname','time new roman')
min_value = max([min(FarSIFs)-0.1, 0]);
max_value = max(FarSIFs)+ 0.1;
plot([min_value max_value], [ b3(1,1) + min_value*b3(2,1) b3(1,1) + max_value*b3(2,1) ], 'r', 'linewidth', 1.5)
axis([0 3 0 3])
set(gca,'fontsize',8,'fontname','time new roman')
box on

% 4
subplot('position',[0.1+0.23*3 0.15 0.18 0.35])
hold on
plot([0 1], [ 0 1 ], 'k', 'linewidth', 0.5);
f1 = scatter(RedSIFs, fitRedSIFs,5,'k','filled')
f1.MarkerFaceAlpha = 0.2;
text(0.05,0.93,['R^2=' num2str(R2s(1,4),'%4.2f')],'fontsize',8,'fontname','time new roman')
text(0.05,0.8,'p < 0.001','fontsize',8,'fontname','time new roman')
min_value = max([min(RedSIFs)-0.1, 0]);
max_value = max(RedSIFs)+ 0.1;
plot([min_value max_value], [ b4(1,1) + min_value*b4(2,1) b4(1,1) + max_value*b4(2,1) ], 'r', 'linewidth', 1.5)
axis([0 1 0 1])
box on
set(gca,'fontsize',8,'fontname','time new roman')

[ax1,h1]=suplabel('Observed Reflectance or SIF (mW m^{-2} nm^{-1} sr^{-1})');
[ax2,h2]=suplabel('Fitted Reflectance or SIF','y');
ax2.Position(1) = ax2.Position(1) - 0.02;
set(h1,'FontSize',10)
set(h2,'FontSize',10)

% save figure
print(gcf, '-dtiff', '-r600', 'figure_S1_fitting.tif')
close all