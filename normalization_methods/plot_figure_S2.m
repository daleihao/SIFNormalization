

clc;
clear all;

%% plot BRDF
figure
set(gcf,'unit','normalized','position',[0.2,0.2,0.5,0.67]);
set(gca, 'Position', [0 0 1 1])
nameStrings = {'Chickpea', 'Grass', 'Rice'};
SZAs = [22 41 21];
for i = 1:3
    load(['results/kk_' nameStrings{i} '.mat']);
    
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
    Relf_NIR_diff = Refl_NIR_est(filters) - Refl_NIR(filters);
    Relf_Red_diff = Refl_Red_est(filters) - Refl_Red(filters) ;
    SIF760_diff = SIF760_est(filters) - SIF760(filters);
    SIF687_diff = SIF687_est(filters) - SIF687(filters);    
      
    if i==1
        istitle  = 1;
    else
        istitle = 0;
    end
    subplot('position',[0.06 0.28*(3-i) + 0.1 0.22 0.26])
        cmin = -max([abs(min(SIF760_diff)) abs(max(SIF760_diff))]);
    cmax = abs(cmin);
cmin = -0.5;
   cmax = abs(cmin);
    plot_polarfigure(RAA_m,VZA_m, SIF760_diff, cmin, cmax, 'Far-red SIF', istitle, SZAs(i))
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
     if( i == 3)
   cb=colorbar('location','southoutside');
    cb.Position = [0.08 0.05 0.18 0.03];
 end 

    subplot('position',[0.29 0.28*(3-i) + 0.1 0.22 0.26])
        cmin = -max([abs(min(SIF687_diff)) abs(max(SIF687_diff))]);
    cmax = abs(cmin);
cmin = -0.401;
   cmax = abs(cmin);

    plot_polarfigure(RAA_m,VZA_m, SIF687_diff, cmin, cmax,'Red SIF', istitle, SZAs(i))
 if( i == 3)
   cb=colorbar('location','southoutside');
    cb.Position = [0.31 0.05 0.18 0.03];
 end 

    
    subplot('position',[0.52 0.28*(3-i) + 0.1 0.22 0.26])    
    cmin = -max([abs(min(Relf_NIR_diff)) abs(max(Relf_NIR_diff))]);
    cmax = abs(cmin); 
    cmin = -0.05;
   cmax = abs(cmin);

    plot_polarfigure(RAA_m,VZA_m, Relf_NIR_diff, cmin, cmax, 'NIR Reflectance', istitle, SZAs(i))
   
 if( i == 3)
   cb=colorbar('location','southoutside');
    cb.Position = [0.54 0.05 0.18 0.03];
 end 

    subplot('position',[0.75 0.28*(3-i) + 0.1 0.22 0.26])  
         cmin = -max([abs(min(Relf_Red_diff)) abs(max(Relf_Red_diff))]);
    cmax = abs(cmin);
    cmin = -0.01;
   cmax = abs(cmin);

    plot_polarfigure(RAA_m,VZA_m, Relf_Red_diff, cmin, cmax, 'Red Reflectance', istitle, SZAs(i))
 if( i == 3)
   cb=colorbar('location','southoutside');
    cb.Position = [0.77 0.05 0.18 0.03];
 end 
end
%%
print(gcf, '-dtiff', '-r600', ['../writting/figure/figureS2.tif'])
close all

%% red
