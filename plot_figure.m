
clc;
clear all;
close all;

%% plot MAE MME
nameStrings = {'Chickpea', 'Grass', 'Rice'};
figure;
for i = 1:3
    load(['results/7_noon_' nameStrings{i} '.mat']);
    DOY_mean = DOY_mean - floor(DOY_mean);
    subplot(2,3,i)
    box on
    hold on
    filters = rMAE760(:,1)<1;
    scatter(DOY_mean(filters), rMAE760(filters,1) *100, 'ko', 'filled')
  %  filters = rMAE760(:,2)<1;
   % scatter(DOY_mean(filters), rMAE760(filters,2) *100, 'go', 'filled')
    filters = rMAE760(:,3)<1;
    scatter(DOY_mean(filters), rMAE760(filters,3) *100, 'bo', 'filled')
    filters = rMAE760(:,4)<1;
    scatter(DOY_mean(filters), rMAE760(filters,4) *100, 'ro', 'filled')
    
    title(nameStrings{i})
    if(i==1)
        ylabel('rMAE (%) for SIF760')
        legend({'Original','NIRv','Kernel-Driven'})
    end
    
    set(gca, 'linewidth', 1, 'xtick', [], 'fontsize', 12)
    xlim(  [min(DOY_mean)- 0.05, max(DOY_mean) + 0.05])
    subplot(2,3,i+3)
    box on
    hold on
    scatter(DOY_mean, rMAE687(:,1) *100, 'ko', 'filled')
    %scatter(DOY_mean, rMAE687(:,2) *100, 'go', 'filled')
    scatter(DOY_mean, rMAE687(:,3) *100, 'bo', 'filled')
    scatter(DOY_mean, rMAE687(:,4) *100, 'ro', 'filled') %R2_687*100,
    xlabel('Hour of Day')
    if(i==1)
        ylabel('rMAE (%) for SIF687')
        legend({'Original','Redv','Kernel-Driven'})
    end
        xlim(  [min(DOY_mean)- 0.05, max(DOY_mean) + 0.05])

    %set(gca)
    set(gca, 'linewidth', 1, 'fontsize', 12)
end

%% plot ANXI
nameStrings = {'Chickpea', 'Grass', 'Rice'};
figure;
for i = 1:3
    load(['results/k_' nameStrings{i} '.mat']);
    subplot(2,3,i)
    box on
    hold on
    filters = ANXI_760(:,1)<100;
    scatter(DOY_mean(filters), ANXI_760(filters,1), 'ko', 'filled')
    filters = ANXI_760(:,2)<100;
    scatter(DOY_mean(filters), ANXI_760(filters,2), 'go', 'filled')
    filters = ANXI_760(:,3)<100;
    scatter(DOY_mean(filters), ANXI_760(filters,3), 'bo', 'filled')
    filters = ANXI_760(:,4)<100;
    scatter(DOY_mean(filters), ANXI_760(filters,4), 'ro', 'filled')
    
    title(nameStrings{i})
    if(i==1)
        ylabel('ANXI for SIF760')
        legend({'Lambertian','NIR','NIRv','Kernel-Driven'})
    end
    
    set(gca, 'linewidth', 1, 'xtick', [], 'fontsize', 12)
    subplot(2,3,i+3)
    box on
    hold on
     filters = ANXI_687(:,1)<100;
    scatter(DOY_mean(filters), ANXI_687(filters,1), 'ko', 'filled')
       filters = ANXI_687(:,2)<100;
    scatter(DOY_mean(filters), ANXI_687(filters,2), 'go', 'filled')
       filters = ANXI_687(:,3)<100;
    scatter(DOY_mean(filters), ANXI_687(filters,3), 'bo', 'filled')
       filters = ANXI_687(:,4)<100;
    scatter(DOY_mean(filters), ANXI_687(filters,4), 'ro', 'filled')
    xlabel('DOY')
    if(i==1)
        ylabel('ANXI for SIF687')
        legend({'Lambertian','Red','Redv','Kernel-Driven'})
    end
end


%% plot CV
nameStrings = {'Chickpea', 'Grass', 'Rice'};
figure;
for i = 1:3
    load(['results/' nameStrings{i} '.mat']);
    subplot(2,3,i)
    box on
    hold on
    filters = CV_760(:,1)<100;
    scatter(DOY_mean(filters), CV_760(filters,1), 'ko', 'filled')
    filters = CV_760(:,2)<100;
    scatter(DOY_mean(filters), CV_760(filters,2), 'go', 'filled')
    filters = CV_760(:,3)<100;
    scatter(DOY_mean(filters), CV_760(filters,3), 'bo', 'filled')
    filters = CV_760(:,4)<100;
    scatter(DOY_mean(filters), CV_760(filters,4), 'ro', 'filled')
    
    title(nameStrings{i})
    if(i==1)
        ylabel('CV for SIF760')
        legend({'Lambertian','NIR','NIRv','Kernel-Driven'})
    end
    
    set(gca, 'linewidth', 1, 'xtick', [], 'fontsize', 12)
    subplot(2,3,i+3)
    box on
    hold on
     filters = CV_687(:,1)<100;
    scatter(DOY_mean(filters), CV_687(filters,1), 'ko', 'filled')
       filters = CV_687(:,2)<100;
    scatter(DOY_mean(filters), CV_687(filters,2), 'go', 'filled')
       filters = CV_687(:,3)<100;
    scatter(DOY_mean(filters), CV_687(filters,3), 'bo', 'filled')
       filters = CV_687(:,4)<100;
    scatter(DOY_mean(filters), CV_687(filters,4), 'ro', 'filled')
    xlabel('DOY')
    if(i==1)
        ylabel('CV for SIF687')
        legend({'Lambertian','Red','Redv','Kernel-Driven'})
    end
end









%% plot
figure;
hold on
plot(nSIF760_true, 'r')
plot(rSIF760, 'g')
plot(nSIF760_NIR, 'b')
plot(nSIF760_NIRv, 'k')
plot(nSIF760_Kernel, 'y')

%% hist
all_differences = [rSIF760-nSIF760_true  nSIF760_NIR-nSIF760_true nSIF760_NIRv-nSIF760_true nSIF760_Kernel-nSIF760_true];
figure;
subplot(1,3,1)
hist(all_differences(VZA<=15, :))
subplot(1,3,2)
hist(all_differences(VZA>15 & VZA<=30, :))
subplot(1,3,3)
hist(all_differences(VZA>30, :))

%%
figure;
subplot(1,2,1)
hold on
scatter(SZA, rSIF760-nSIF760_true, 'r')
scatter(SZA, nSIF760_NIR-nSIF760_true, 'g')
scatter(SZA, nSIF760_NIRv-nSIF760_true, 'b')
scatter(SZA, nSIF760_Kernel-nSIF760_true, 'y')


subplot(1,2,2)
hold on
scatter(VZA, rSIF687-nSIF687_true, 'r')
scatter(VZA, nSIF687_Red-nSIF687_true, 'g')
scatter(VZA, nSIF687_Redv-nSIF687_true, 'b')
scatter(VZA, nSIF687_Kernel-nSIF687_true, 'y')

disp(corrcoef(rSIF760, nSIF760_true))
disp(corrcoef(nSIF760_NIR, nSIF760_true))
disp(corrcoef(nSIF760_NIRv, nSIF760_true))
disp(corrcoef(nSIF760_Kernel, nSIF760_true))

disp(corrcoef(rSIF687, nSIF687_true))
disp(corrcoef(nSIF687_Red, nSIF687_true))
disp(corrcoef(nSIF687_Redv, nSIF687_true))
disp(corrcoef(nSIF687_Kernel, nSIF687_true))