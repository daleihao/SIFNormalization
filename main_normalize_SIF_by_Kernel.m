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
    rSIF687 = SIF687%./PAR*1e2;
    rSIF760 = SIF760%./PAR*1e2;
    
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
    
    Refl_NIR_est = zeros(counts, 1);
    Refl_Red_est = zeros(counts, 1);
    SIF760_est = zeros(counts, 1);
    SIF687_est = zeros(counts, 1);
    
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
    R2_NIR = zeros(cycle_n, 1);
    R2_Red = zeros(cycle_n, 1);
    RMSE_760 = zeros(cycle_n, 1);
    RMSE_687 = zeros(cycle_n, 1);
    RMSE_NIR = zeros(cycle_n, 1);
    RMSE_Red = zeros(cycle_n, 1);
    rRMSE_760 = zeros(cycle_n, 1);
    rRMSE_687 = zeros(cycle_n, 1);
    rRMSE_NIR = zeros(cycle_n, 1);
    rRMSE_Red = zeros(cycle_n, 1);

    DOY_mean = zeros(cycle_n, 1);
    
    ANXI_760 = zeros(cycle_n, 4);
    ANXI_687 = zeros(cycle_n, 4);
    CV_760 = zeros(cycle_n, 4);
    CV_687 = zeros(cycle_n, 4);
    for cycle_i = 1:cycle_n
        
        filters = Cycle_Num == cycle_i  & SIF687>0.2 & SIF760>0.2;
        filters_nadir = Cycle_Num == cycle_i & VZA == 0 & SIF687>0.2 & SIF760>0.2;
        count_n = sum(filters);
        
        [Kparms_NIR, R2, RMSE, Relf_NIR_est_tmp] = kernelParameterRetrieval_refl(SZA(filters), VZA(filters), RAA(filters), Refl_NIR(filters), [nameStrings{name_i} '_NIR_' num2str(cycle_i)], 0);
        R2_NIR(cycle_i) = R2;
        RMSE_NIR(cycle_i) = RMSE;
        [Kparms_Red, R2, RMSE, Relf_Red_est_tmp] = kernelParameterRetrieval_refl(SZA(filters), VZA(filters), RAA(filters), Refl_Red(filters), [nameStrings{name_i} '_Red_' num2str(cycle_i)], 0);
        R2_Red(cycle_i) = R2;
        RMSE_Red(cycle_i) = RMSE;
        

        kNIR  = CalculateRefl(SZA(filters), VZA(filters), RAA(filters), Kparms_NIR);
        kRed  = CalculateRefl(SZA(filters), VZA(filters), RAA(filters), Kparms_Red);
        kNIR0  = CalculateRefl(SZA(filters), zeros(count_n, 1), zeros(count_n, 1), Kparms_NIR);
        kRed0  = CalculateRefl(SZA(filters), zeros(count_n, 1), zeros(count_n, 1), Kparms_Red);
       
        kNDVI = (kNIR-kRed)./(kNIR+kRed);
        kNDVI0 = (kNIR0-kRed0)./(kNIR0+kRed0);
        kNIRv = kNDVI.*kNIR;
        kNIRv0 = kNDVI0.*kNIR0;
        kRedv = kRed.*kNDVI.*kNDVI;
        kRedv0 = kRed0.*kNDVI0.*kNDVI0;
        
        
        
        nSIF760_true(filters) = mean(rSIF760(filters_nadir));
        nSIF760_NIR(filters) = kNIR0./kNIR.*rSIF760(filters);
        nSIF760_NIRv(filters) = kNIRv0./kNIRv.*rSIF760(filters);
        
        nSIF687_true(filters) = mean(rSIF687(filters_nadir));
        nSIF687_Red(filters) = kRed0./kRed.*rSIF687(filters);
        nSIF687_Redv(filters) = kRedv0./kRedv.*rSIF687(filters);
        %nSIF687_EVI2(filters) = mean(EVI2(filters_nadir))./EVI2(filters).*rSIF687(filters);
        
        %% normalized by kernel-driven models
        [Kparms, R2, RMSE, SIF760_est_tmp] = kernelParameterRetrieval_refl(SZA(filters), VZA(filters), RAA(filters), rSIF760(filters), [nameStrings{name_i} '_760_' num2str(cycle_i)], 0);
        c_factor = CalculateCorrectionFactor(SZA(filters), VZA(filters), RAA(filters), Kparms);
        nSIF760_Kernel(filters) = c_factor.*rSIF760(filters);
        
        R2_760(cycle_i) = R2;
        RMSE_760(cycle_i) = RMSE;
        
        [Kparms, R2, RMSE, SIF687_est_tmp] = kernelParameterRetrieval_refl(SZA(filters), VZA(filters), RAA(filters), rSIF687(filters), [nameStrings{name_i} '_687_' num2str(cycle_i)], 0);
        c_factor = CalculateCorrectionFactor(SZA(filters), VZA(filters), RAA(filters), Kparms);
        nSIF687_Kernel(filters) = c_factor.*rSIF687(filters);
        
        R2_687(cycle_i) = R2;
        RMSE_687(cycle_i) = RMSE;
        
        
    Refl_NIR_est(filters) = Relf_NIR_est_tmp;
    Refl_Red_est(filters) = Relf_Red_est_tmp;
    SIF760_est(filters) = SIF760_est_tmp;
    SIF687_est(filters) = SIF687_est_tmp;
        %% accuracy statistics
        SZA_mean(cycle_i, 1) = mean(SZA(filters));
        DOY_mean(cycle_i, 1) = mean(DOY(filters));
        
    rRMSE_760(cycle_i, 1) = RMSE_760(cycle_i, 1)/nanmean(rSIF760(filters));
    rRMSE_687(cycle_i, 1) = RMSE_687(cycle_i, 1)/nanmean(rSIF687(filters));
    rRMSE_NIR(cycle_i, 1) = RMSE_NIR(cycle_i, 1)/nanmean(Refl_NIR(filters));
    rRMSE_Red(cycle_i, 1) = RMSE_Red(cycle_i, 1)/nanmean(Refl_Red(filters));

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
    save(['results/kk_' nameStrings{name_i} '.mat']);
end