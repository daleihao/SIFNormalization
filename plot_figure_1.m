close all
figure;
set(gcf,'unit','normalized','position',[0.2,0.2,0.7,0.6]);
set(gca, 'Position', [0 0 1 1])

nameStrings = {'Chickpea', 'Grass', 'Rice'};
for i = 1:3
    load(['results/k_' nameStrings{i} '.mat']);
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
    subplot('position',[0.065 + (0.25 + 0.06)*(i-1) 0.55 0.25 0.35])
    hold on
    
    plot(DOY_mean, R2_760, 'bo-', 'MarkerFaceColor', 'b')
    plot(DOY_mean, R2_687, 'b^-', 'MarkerFaceColor', 'b')
    plot(DOY_mean, R2_NIR, 'ro-', 'MarkerFaceColor', 'r')
    plot(DOY_mean, R2_Red, 'r^-', 'MarkerFaceColor', 'r')
    if(i == 1)
        ylabel('R^2')
    end
    axis([ min(DOY_mean)-0.5 max(DOY_mean)+0.5 0 1])
    box on
    set(gca, 'linewidth', 1, 'fontsize', 12) % 'xticklabel', {}
    if(i == 1)
        legend({'Far-red SIF', 'Red SIF', 'NIR Reflectance', 'Red Reflectance'}, 'Location','southeast');
    end
    title(nameStrings{i}, 'FontWeight', 'Bold')
    subplot('position',[0.065 + (0.25 + 0.06)*(i-1) 0.12 0.25 0.35])
    hold on
    
    [hAx_1,hLine1_1,hLine2_1] = plotyy(DOY_mean, rRMSE_760, DOY_mean, rRMSE_NIR);
    [hAx_2,hLine1_2,hLine2_2] =  plotyy(DOY_mean, rRMSE_687, DOY_mean, rRMSE_Red);
    
    set(hAx_1,{'ycolor'},{'b';'r'})
    set(hAx_1,{'linewidth'},{1; 1})
    set(hAx_2,{'linewidth'},{1; 1})
    
    set(hAx_2,{'ycolor'},{'b';'r'})
    set(hAx_1, {'xlim'}, {   [min(DOY_mean)- 0.5, max(DOY_mean) + 0.5];[min(DOY_mean)- 0.5, max(DOY_mean)+ 0.5]});
    set(hAx_2, {'xlim'}, {   [min(DOY_mean)- 0.5, max(DOY_mean)+ 0.5] ;[min(DOY_mean)- 0.5, max(DOY_mean)+ 0.5]});
    set(hAx_1, {'ylim'}, {   [0 0.45];[0 0.15]});
    set(hAx_2, {'ylim'}, {   [0 0.45];[0 0.15]});
    set(hAx_1, {'ytick'}, {   [0 0.2 0.4];[0 0.05 0.1 0.15]});
    set(hAx_2, {'ytick'}, {   [0 0.2 0.4];[0 0.05 0.1 0.15]});
    set(hAx_1, {'yTickLabel'}, {   {};{}});
    set(hAx_2, {'yTickLabel'}, {   {'0', '20', '40'};{'0', '5' ,'10', '15'}});
    
    set(hAx_1, {'fontSize'}, {   12;12});
    set(hAx_2, {'fontSize'}, {   12;12});
    
    
    hLine1_1.MarkerFaceColor = 'b';
    hLine1_1.MarkerEdgeColor = 'b';
    hLine1_1.Marker = 'o';
    hLine1_1.LineStyle = '-';
    hLine1_1.Color = 'b';
     
    hLine2_1.MarkerFaceColor = 'r';
    hLine2_1.MarkerEdgeColor = 'r';
    hLine2_1.Marker = 'o';
    hLine2_1.LineStyle = '-';
        hLine2_1.Color = 'r';

    hLine1_2.MarkerFaceColor = 'b';
    hLine1_2.MarkerEdgeColor = 'b';
    hLine1_2.Marker = '^';
    hLine1_2.LineStyle = '-';
        hLine1_2.Color = 'b';

    hLine2_2.MarkerFaceColor = 'r';
    hLine2_2.MarkerEdgeColor = 'r';
    hLine2_2.Marker = '^';
    hLine2_2.LineStyle = '-';
        hLine2_2.Color = 'r';

    box on
    if i == 1
        ylabel( hAx_1(1), sprintf('%s','rRMSE for SIF (%)'), 'fontsize', 12) %r\n10^{-2}m^{-2}um^{-1}
    end
    if i == 3
        ylabel( hAx_2(2), 'rRMSE for Reflectance (%)', 'fontsize', 12);
    end
    
    xlabel('Hour of Day')
    set(gca, 'fontsize', 12)
end
print(gcf, '-dtiff', '-r600', ['../writting/figure/figure1.tif'])
close all
