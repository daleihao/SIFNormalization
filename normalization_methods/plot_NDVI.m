

clc;
clear all;

%% plot BRDF
figure
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
    NDVI_m = NDVI(filters);
    NIRv_m = NIRv(filters);
    NIR_m = Refl_NIR(filters);
    Red_m = Refl_Red(filters);

    SIF_m = rSIF760(filters);
    cmin = 0;
    cmax = 1;
    subplot(3,5,1+(i-1)*5)
    cmin = min(NDVI_m);
    cmax = max(NDVI_m);
    plot_polarfigure(VAA_m,VZA_m, NDVI_m, cmin, cmax, 'NDVI')
    subplot(3,5,2+(i-1)*5)
    cmin = min(NIRv_m);
    cmax = max(NIRv_m);

    plot_polarfigure(VAA_m,VZA_m, NIRv_m, cmin, cmax, 'NIRv')
    cmin = min(NIR_m);
    cmax = max(NIR_m);
    
    subplot(3,5,3+(i-1)*5)
    plot_polarfigure(VAA_m,VZA_m, NIR_m, cmin, cmax,'NIR')

    cmin = min(SIF_m);
    cmax = max(SIF_m);
    
    subplot(3,5,4+(i-1)*5)
    plot_polarfigure(VAA_m,VZA_m, SIF_m, cmin, cmax,'SIF')

    
    cmin = min(Red_m);
    cmax = max(Red_m);
    
    subplot(3,5,5+(i-1)*5)
    plot_polarfigure(VAA_m,VZA_m, Red_m, cmin, cmax,'Red')

end

