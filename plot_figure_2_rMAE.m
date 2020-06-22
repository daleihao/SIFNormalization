clc;
clear all;
close all;

%% plot MAE MME
nameStrings = {'Chickpea', 'Grass', 'Rice'};
figure;
set(gcf,'unit','normalized','position',[0.2,0.2,0.7,0.6]);
set(gca, 'Position', [0 0 1 1])

for i = 1:3
    load(['results/no_nadir_' nameStrings{i} '.mat']);
    DOY_mean = DOY_mean - floor(DOY_mean);    
    
      % from UTC to local
    if(i == 1)
        DOY_mean = DOY_mean + 11.07836/15/24;
    elseif (i == 2)
        DOY_mean = DOY_mean - 5.77913/15/24;
    else
        DOY_mean = DOY_mean + 11.06905/15/24;
    end;
    DOY_mean = DOY_mean*24;
    %subplot('position',[0.03 + 0.32*(i-1) 0.55 0.3 0.35])
    subplot('position',[0.07 + (0.28 + 0.04)*(i-1) 0.55 0.28 0.35])
    
    
    box on
    hold on
    filters = rMAE760(:,1)<1;
    plot(DOY_mean(filters), rMAE760(filters,1) *100, 'ko-', 'MarkerFaceColor', 'k')
  %  filters = rMAE760(:,2)<1;
   % scatter(DOY_mean(filters), rMAE760(filters,2) *100, 'go', 'filled')
    filters = rMAE760(:,3)<1;
    plot(DOY_mean(filters), rMAE760(filters,3) *100, 'bo-', 'MarkerFaceColor', 'b')
    filters = rMAE760(:,4)<1;
    plot(DOY_mean(filters), rMAE760(filters,4) *100, 'ro-', 'MarkerFaceColor', 'r')
    
    title(nameStrings{i}, 'FontWeight', 'Bold')
    if(i==1)
        ylabel('rMAE (%) for far-red SIF')
        legend({'Raw Data','NIRv/Redv-based','KD-based'})
    end
    
    set(gca, 'linewidth', 1, 'fontsize', 12) % 'xtick', []
    xlim(  [min(DOY_mean)- 0.5, max(DOY_mean) + 0.5])
     ylim( [0 70])
    subplot('position',[0.07 + (0.28 + 0.04)*(i-1) 0.12 0.28 0.35])
    box on
    hold on
    filters = rMAE687(:,1)>=0;
    plot(DOY_mean(filters), rMAE687(filters,1) *100, 'ko-', 'MarkerFaceColor', 'k')
    %scatter(DOY_mean, rMAE687(:,2) *100, 'go', 'filled')
    plot(DOY_mean(filters), rMAE687(filters,3) *100, 'bo-', 'MarkerFaceColor', 'b')
    plot(DOY_mean(filters), rMAE687(filters,4) *100, 'ro-', 'MarkerFaceColor', 'r') %R2_687*100,
    xlabel('Hour of Day')
    if(i==1)
        ylabel('rMAE (%) for red SIF')
       % legend({'Original','Redv','Kernel-Driven'})
    end
        xlim(  [min(DOY_mean)- 0.5, max(DOY_mean) + 0.5])
 ylim([ 0 100])
    %set(gca)
    set(gca, 'linewidth', 1, 'fontsize', 12)
end
print(gcf, '-dtiff', '-r600', ['../writting/figure/figure2c_no_nadir.tif'])
close all
