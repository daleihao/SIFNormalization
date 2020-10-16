%%% normalize SIF using different methods
%%% author: Dalei Hao
%%% 2020/4/27

clc;
clear all;
close all;

%% import data
nameStrings = {'Chickpea', 'Grass', 'Rice'};
for name_i = 1:3
    filename = [nameStrings{name_i} '_multi_angular_data.csv'];
    
    allData = importdata(['../data/formatting/' filename], ',', 1);
    textData = allData.textdata;
    allData = allData.data;
    
    Cycle_Num = allData(:,1);
    DOY =  allData(:,2);
    
    VZA = allData(:,3);
    RAA = allData(:,4);
    SZA = allData(:,5);
    SAA = allData(:,6);
    
    PAR = allData(:,7);
    Refl_Red = allData(:,8);
    Refl_NIR = allData(:,9);
    SIF760 = allData(:,10);
    SIF687 = allData(:,11);
    
    
    %% calculate corresponding NIRv, NDVI, EVI2, NIR, rSIF
    rSIF687 = SIF687;%./PAR*1e3;
    rSIF760 = SIF760;%./PAR*1e3;
    
    NDVI = (Refl_NIR-Refl_Red)./(Refl_NIR+Refl_Red);
    %EVI = 2.5*(Refl_NIR-Refl_Red)./(1+Refl_NIR+6*Refl_Red-7.5*Refl_Blue);
    EVI2 =  2.5*(Refl_NIR-Refl_Red)./(1+Refl_NIR+2.4*Refl_Red);
    NIRv = Refl_NIR.*NDVI;
    Redv = Refl_Red.*NDVI.*NDVI;
    
    %% Red
    % SIF normalization by Reflectance
    counts = size(NIRv, 1);
    nSIF760_true = zeros(counts, 1);
    nSIF760_NIR = zeros(counts, 1);
    nSIF760_NIRv = zeros(counts, 1);
    nSIF760_Kernel = zeros(counts, 1);
    nSIF687_true = zeros(counts, 1);
    nSIF687_Red = zeros(counts, 1);
    nSIF687_Redv = zeros(counts, 1);
    nSIF687_Kernel = zeros(counts, 1);
    
    cycle_n = max(Cycle_Num);
    MAE760 = zeros(cycle_n, 4);
    MAE687 = zeros(cycle_n, 4);
    rMAE760 = zeros(cycle_n, 4);
    rMAE687 = zeros(cycle_n, 4);
    MME760 = zeros(cycle_n, 4);
    MME687 = zeros(cycle_n, 4);
    rMME760 = zeros(cycle_n, 4);
    rMME687 = zeros(cycle_n, 4);
    SZA_mean = zeros(cycle_n, 1);
    R2_760 = zeros(cycle_n, 1);
    R2_687 = zeros(cycle_n, 1);
    RMSE_760 = zeros(cycle_n, 1);
    RMSE_687 = zeros(cycle_n, 1);

    DOY_mean = zeros(cycle_n, 1);
    
    ANXI_760 = zeros(cycle_n, 4);
    ANXI_687 = zeros(cycle_n, 4);

    CV_760 = zeros(cycle_n, 4);
    CV_687 = zeros(cycle_n, 4);
    
    R2_NIRv = zeros(cycle_n, 1);
    R2_NIR = zeros(cycle_n, 1);

    for cycle_i = 1:cycle_n
        
        filters = Cycle_Num == cycle_i & rSIF687>0.2 & rSIF760>0.2;
        filters_nadir = Cycle_Num == cycle_i & VZA == 0  & rSIF687>0.2 & rSIF760>0.2;
        
        nSIF760_true(filters) = mean(rSIF760(filters_nadir));
        nSIF760_NIR(filters) = mean(Refl_NIR(filters_nadir))./Refl_NIR(filters).*rSIF760(filters);
        nSIF760_NIRv(filters) = mean(NIRv(filters_nadir))./NIRv(filters).*rSIF760(filters);
        %nSIF760_EVI2(filters) = mean(EVI2(filters_nadir))./EVI2(filters).*rSIF760(filters);
        
        nSIF687_true(filters) = mean(rSIF687(filters_nadir));
        nSIF687_Red(filters) = mean(Refl_Red(filters_nadir))./Refl_Red(filters).*rSIF687(filters);
        nSIF687_Redv(filters) = mean(Redv(filters_nadir))./Redv(filters).*rSIF687(filters);
        %nSIF687_EVI2(filters) = mean(EVI2(filters_nadir))./EVI2(filters).*rSIF687(filters);
        
        figure;
        subplot(1,2,1)
        scatter(NIRv(filters),rSIF760(filters)); 
        R2 = corrcoef(NIRv(filters),rSIF760(filters));
        R2 = R2(1,2);
        R2 = R2^2;
        R2_NIRv(cycle_i, 1) = R2;

        title(['NIRv-' num2str(R2,'%.2f')])
        subplot(1,2,2)
        scatter(Refl_NIR(filters),rSIF760(filters)); 
        title('NIR')
        R2 = corrcoef(Refl_NIR(filters),rSIF760(filters));
        R2 = R2(1,2);
        R2 = R2^2;
        R2_NIR(cycle_i, 1) = R2;
        title(['NIR-' num2str(R2,'%.2f')])
       % print(gcf, '-dtiff', '-r600',['results/NIRv_sif/' nameStrings{name_i} '_' num2str(cycle_i) '.tif']); 
       % close all;
    end
    save(['results/NIRv_sif/' nameStrings{name_i}  '.mat'], 'R2_NIR', 'R2_NIRv')
   
end