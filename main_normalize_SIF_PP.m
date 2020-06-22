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
    
    % select CPP
    RAA = allData(:,4);
    VZA = allData(:,3);
    allData = allData((RAA==90 | RAA==270)|(VZA==0), :);
    
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

    for cycle_i = 1:cycle_n
        
        filters = Cycle_Num == cycle_i & rSIF687>0.2 & rSIF760>0.2;
        filters_nadir = Cycle_Num == cycle_i & VZA == 0  & rSIF687>0.2 & rSIF760>0.2;
        if sum(filters)<=1
        continue;
        end
        nSIF760_true(filters) = mean(rSIF760(filters_nadir));
        nSIF760_NIR(filters) = mean(Refl_NIR(filters_nadir))./Refl_NIR(filters).*rSIF760(filters);
        nSIF760_NIRv(filters) = mean(NIRv(filters_nadir))./NIRv(filters).*rSIF760(filters);
        %nSIF760_EVI2(filters) = mean(EVI2(filters_nadir))./EVI2(filters).*rSIF760(filters);
        
        nSIF687_true(filters) = mean(rSIF687(filters_nadir));
        nSIF687_Red(filters) = mean(Refl_Red(filters_nadir))./Refl_Red(filters).*rSIF687(filters);
        nSIF687_Redv(filters) = mean(Redv(filters_nadir))./Redv(filters).*rSIF687(filters);
        %nSIF687_EVI2(filters) = mean(EVI2(filters_nadir))./EVI2(filters).*rSIF687(filters);
        
        %% normalized by kernel-driven models
        [Kparms, R2, RMSE] = kernelParameterRetrieval(SZA(filters), VZA(filters), RAA(filters), rSIF760(filters),[nameStrings{name_i} '_760_' num2str(cycle_i)], 0);
        c_factor = CalculateCorrectionFactor(SZA(filters), VZA(filters), RAA(filters), Kparms);
        nSIF760_Kernel(filters) = c_factor.*rSIF760(filters);
        
        R2_760(cycle_i) = R2;
        RMSE_760(cycle_i) = RMSE;

        
        [Kparms, R2, RMSE] = kernelParameterRetrieval(SZA(filters), VZA(filters), RAA(filters), rSIF687(filters), [nameStrings{name_i} '_687_' num2str(cycle_i)], 0);
        c_factor = CalculateCorrectionFactor(SZA(filters), VZA(filters), RAA(filters), Kparms);
        nSIF687_Kernel(filters) = c_factor.*rSIF687(filters);
        
        R2_687(cycle_i) = R2;
        RMSE_687(cycle_i) = RMSE;

        %% accuracy statistics
        SZA_mean(cycle_i, 1) = mean(SZA(filters));
        DOY_mean(cycle_i, 1) = mean(DOY(filters));
        
        % NIR
        MAE760(cycle_i, 1) = mean(abs(nSIF760_true(filters)- rSIF760(filters)));
        MAE760(cycle_i, 2) = mean(abs(nSIF760_true(filters)- nSIF760_NIR(filters)));
        MAE760(cycle_i, 3) = mean(abs(nSIF760_true(filters)- nSIF760_NIRv(filters)));
        MAE760(cycle_i, 4) = mean(abs(nSIF760_true(filters)- nSIF760_Kernel(filters)));
        
        MME760(cycle_i, 1) = max(abs(nSIF760_true(filters)- rSIF760(filters)));
        MME760(cycle_i, 2) = max(abs(nSIF760_true(filters)- nSIF760_NIR(filters)));
        MME760(cycle_i, 3) = max(abs(nSIF760_true(filters)- nSIF760_NIRv(filters)));
        MME760(cycle_i, 4) = max(abs(nSIF760_true(filters)- nSIF760_Kernel(filters)));
        
        rMAE760(cycle_i, 1) = mean(abs(nSIF760_true(filters)- rSIF760(filters)))/mean(nSIF760_true(filters));
        rMAE760(cycle_i, 2) = mean(abs(nSIF760_true(filters)- nSIF760_NIR(filters)))/mean(nSIF760_true(filters));
        rMAE760(cycle_i, 3) = mean(abs(nSIF760_true(filters)- nSIF760_NIRv(filters)))/mean(nSIF760_true(filters));
        rMAE760(cycle_i, 4) = mean(abs(nSIF760_true(filters)- nSIF760_Kernel(filters)))/mean(nSIF760_true(filters));
        
        rMME760(cycle_i, 1) = max(abs(nSIF760_true(filters)- rSIF760(filters)))/mean(nSIF760_true(filters));
        rMME760(cycle_i, 2) = max(abs(nSIF760_true(filters)- nSIF760_NIR(filters)))/mean(nSIF760_true(filters));
        rMME760(cycle_i, 3) = max(abs(nSIF760_true(filters)- nSIF760_NIRv(filters)))/mean(nSIF760_true(filters));
        rMME760(cycle_i, 4) = max(abs(nSIF760_true(filters)- nSIF760_Kernel(filters)))/mean(nSIF760_true(filters));
        
        % red
        MAE687(cycle_i, 1) = mean(abs(nSIF687_true(filters)- rSIF687(filters)));
        MAE687(cycle_i, 2) = mean(abs(nSIF687_true(filters)- nSIF687_Red(filters)));
        MAE687(cycle_i, 3) = mean(abs(nSIF687_true(filters)- nSIF687_Redv(filters)));
        MAE687(cycle_i, 4) = mean(abs(nSIF687_true(filters)- nSIF687_Kernel(filters)));
        
        MME687(cycle_i, 1) = max(abs(nSIF687_true(filters)- rSIF687(filters)));
        MME687(cycle_i, 2) = max(abs(nSIF687_true(filters)- nSIF687_Red(filters)));
        MME687(cycle_i, 3) = max(abs(nSIF687_true(filters)- nSIF687_Redv(filters)));
        MME687(cycle_i, 4) = max(abs(nSIF687_true(filters)- nSIF687_Kernel(filters)));
        
        rMAE687(cycle_i, 1) = mean(abs(nSIF687_true(filters)- rSIF687(filters)))/mean(nSIF687_true(filters));
        rMAE687(cycle_i, 2) = mean(abs(nSIF687_true(filters)- nSIF687_Red(filters)))/mean(nSIF687_true(filters));
        rMAE687(cycle_i, 3) = mean(abs(nSIF687_true(filters)- nSIF687_Redv(filters)))/mean(nSIF687_true(filters));
        rMAE687(cycle_i, 4) = mean(abs(nSIF687_true(filters)- nSIF687_Kernel(filters)))/mean(nSIF687_true(filters));
        
        rMME687(cycle_i, 1) = max(abs(nSIF687_true(filters)- rSIF687(filters)))/mean(nSIF687_true(filters));
        rMME687(cycle_i, 2) = max(abs(nSIF687_true(filters)- nSIF687_Red(filters)))/mean(nSIF687_true(filters));
        rMME687(cycle_i, 3) = max(abs(nSIF687_true(filters)- nSIF687_Redv(filters)))/mean(nSIF687_true(filters));
        rMME687(cycle_i, 4) = max(abs(nSIF687_true(filters)- nSIF687_Kernel(filters)))/mean(nSIF687_true(filters));
        
        ANXI_760(cycle_i, 1) = max(rSIF760(filters))/min(rSIF760(filters));
        ANXI_760(cycle_i, 2) = max(nSIF760_NIR(filters))/min(nSIF760_NIR(filters));
        ANXI_760(cycle_i, 3) = max(nSIF760_NIRv(filters))/min(nSIF760_NIRv(filters));
        ANXI_760(cycle_i, 4) = max(nSIF760_Kernel(filters))/min(nSIF760_Kernel(filters));
        
        ANXI_687(cycle_i, 1) = max(rSIF687(filters))/min(rSIF687(filters));
        ANXI_687(cycle_i, 2) = max(nSIF687_Red(filters))/min(nSIF687_Red(filters));
        ANXI_687(cycle_i, 3) = max(nSIF687_Redv(filters))/min(nSIF687_Redv(filters));
        ANXI_687(cycle_i, 4) = max(nSIF687_Kernel(filters))/min(nSIF687_Kernel(filters));
  
        CV_760(cycle_i, 1) = std(rSIF760(filters))/mean(rSIF760(filters));
        CV_760(cycle_i, 2) = std(nSIF760_NIR(filters))/mean(nSIF760_NIR(filters));
        CV_760(cycle_i, 3) = std(nSIF760_NIRv(filters))/mean(nSIF760_NIRv(filters));
        CV_760(cycle_i, 4) = std(nSIF760_Kernel(filters))/mean(nSIF760_Kernel(filters));
        
        CV_687(cycle_i, 1) = std(rSIF687(filters))/mean(rSIF687(filters));
        CV_687(cycle_i, 2) = std(nSIF687_Red(filters))/mean(nSIF687_Red(filters));
        CV_687(cycle_i, 3) = std(nSIF687_Redv(filters))/mean(nSIF687_Redv(filters));
        CV_687(cycle_i, 4) = std(nSIF687_Kernel(filters))/mean(nSIF687_Kernel(filters));

    end
   
   %% save 
    save(['results/CPP_' nameStrings{name_i} '.mat']);
end