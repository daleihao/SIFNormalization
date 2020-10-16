

clc;
clear all;
SZAs = [22 41 21];
%% plot BRDF
figure
set(gcf,'unit','normalized','position',[0.2,0.2,0.38,0.64]);
set(gca, 'Position', [0 0 1 1])
  htitle =  suptitle('Far-red SIF');
  htitle.FontWeight = 'Bold';

nameStrings = {'Chickpea', 'Grass', 'Rice'};
for i = 1:3
    load(['results/' nameStrings{i} '.mat']);
    
    switch i
        case 1
            index = 7; %8 7 and 8
        case 2
            index = 11; % 9 and 11
        case 3
            index = 7;% 4, 7, 9, 13
    end
    filters = Cycle_Num == index & rSIF760>0.2 & rSIF687>0.2;
    RAA_m = RAA(filters);
    VAA_m = SAA(filters)+RAA(filters);
    VAA_m(VAA_m>360)=VAA_m(VAA_m>360)-360;
    VZA_m = VZA(filters);
    rSIF760_m = rSIF760(filters) - nSIF760_true(filters);
    nSIF760_NIR_m = nSIF760_NIR(filters) - nSIF760_true(filters);
    nSIF760_NIRv_m = nSIF760_NIRv(filters) - nSIF760_true(filters);
    nSIF760_Kernel_m = nSIF760_Kernel(filters) - nSIF760_true(filters);
    
    cmin = -max([abs(min(rSIF760_m)) abs(max(rSIF760_m))]);
    cmax = abs(cmin);
    cmin = -0.7;
    cmax = 0.7;
    
    if i==1
        istitle  = 1;
    else
        istitle = 0;
    end
    subplot('position',[0.06 0.29*(3-i) + 0.01 0.25 0.27])
    plot_polarfigure(RAA_m,VZA_m, rSIF760_m, cmin, cmax, 'Raw Data', istitle, SZAs(i))
    % ylabel('Chickpea')
       if(i == 1)
        x1h = ylabel('Chickpea', 'FontWeight', 'Bold')
        x1h.Position(1) =   x1h.Position(1) + 15; 
    elseif(i == 2)
        x1h = ylabel('Grass', 'FontWeight', 'Bold')
        x1h.Position(1) =   x1h.Position(1) + 15; 
    else
        x1h = ylabel('Rice', 'FontWeight', 'Bold')
        x1h.Position(1) =   x1h.Position(1) + 15; 
    end
    subplot('position',[0.35 0.29*(3-i) + 0.01 0.25 0.27])
    plot_polarfigure(RAA_m,VZA_m, nSIF760_NIRv_m, cmin, cmax,'NIRv-based', istitle, SZAs(i))
    
    subplot('position',[0.64 0.29*(3-i) + 0.01 0.25 0.27])
    
    plot_polarfigure(RAA_m,VZA_m, nSIF760_Kernel_m, cmin, cmax, 'KD-based', istitle, SZAs(i))
    cb=colorbar('location','eastoutside');
    cb.Position = [0.91 0.29*(3-i)+0.03 0.02 0.25];
end

print(gcf, '-dtiff', '-r600', ['../writting/figure/figure2a.tif'])
close all

%% red
%% plot BRDF
figure
set(gcf,'unit','normalized','position',[0.2,0.2,0.38,0.64]);
set(gca, 'Position', [0 0 1 1])

  htitle =  suptitle('Red SIF');
  htitle.FontWeight = 'Bold';

nameStrings = {'Chickpea', 'Grass', 'Rice'};
for i = 1:3
    load(['results/' nameStrings{i} '.mat']);
    
    switch i
        case 1
            index = 7; %8 7 and 8
        case 2
            index = 11; % 9 and 11
        case 3
            index = 7;% 4, 7, 9, 13
    end
    filters = Cycle_Num == index & rSIF687>0.2 & rSIF687>0.2;
    RAA_m = RAA(filters);
    VAA_m = SAA(filters)+RAA(filters);
    VAA_m(VAA_m>360)=VAA_m(VAA_m>360)-360;
    VZA_m = VZA(filters);
    rSIF687_m = rSIF687(filters) - nSIF687_true(filters);
    nSIF687_Red_m = nSIF687_Red(filters) - nSIF687_true(filters);
    nSIF687_Redv_m = nSIF687_Redv(filters) - nSIF687_true(filters);
    nSIF687_Kernel_m = nSIF687_Kernel(filters) - nSIF687_true(filters);
    
    cmin = -max([abs(min(rSIF687_m)) abs(max(rSIF687_m))]);
    cmax = abs(cmin);
    
    cmin = -0.7;
    cmax = 0.7;

    if i==1
        istitle  = 1;
    else
        istitle = 0;
    end
    subplot('position',[0.06 0.29*(3-i) + 0.01 0.25 0.27])
    plot_polarfigure(RAA_m,VZA_m, rSIF687_m, cmin, cmax, 'Raw Data', istitle, SZAs(i))
    % ylabel('Chickpea')
       if(i == 1)
        x1h = ylabel('Chickpea', 'FontWeight', 'Bold')
        x1h.Position(1) =   x1h.Position(1) + 15; 
    elseif(i == 2)
        x1h = ylabel('Grass', 'FontWeight', 'Bold')
        x1h.Position(1) =   x1h.Position(1) + 15; 
    else
        x1h = ylabel('Rice', 'FontWeight', 'Bold')
        x1h.Position(1) =   x1h.Position(1) + 15; 
    end
    subplot('position',[0.35 0.29*(3-i) + 0.01 0.25 0.27])
    plot_polarfigure(RAA_m,VZA_m, nSIF687_Redv_m, cmin, cmax,'Redv-based', istitle, SZAs(i))
    
    subplot('position',[0.64 0.29*(3-i) + 0.01 0.25 0.27])
    
    plot_polarfigure(RAA_m,VZA_m, nSIF687_Kernel_m, cmin, cmax, 'KD-based', istitle, SZAs(i))
    cb=colorbar('location','eastoutside');
    cb.Position = [0.91 0.29*(3-i)+0.03 0.02 0.25];
            
end

print(gcf, '-dtiff', '-r600', ['../writting/figure/figure2b.tif'])
close all

