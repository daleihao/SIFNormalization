%%% normalize SIF using different methods
%%% author: Dalei Hao
%%% 2020/4/27

clc;
clear all;
close all;

%% import data
nameStrings = {'Chickpea', 'Grass', 'Rice'};
for obs_numbers = 4:16
    
    for name_i = 1:3
       
            rMAE760 = zeros(100, 4);
            rMAE687 = zeros(100, 4);
        for rand_n = 1:100
            disp(rand_n)

            filename = [nameStrings{name_i} '_multi_angular_data.csv'];
            
            allData = importdata(['../data/formatting/' filename], ',', 1);
            textData = allData.textdata;
            allData = allData.data;
            
            SIF760 = allData(:,10);
            SIF687 = allData(:,11);
            filters = SIF687>0.2 & SIF760>0.2;
            allData = allData(filters, :);
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
            rSIF687 = SIF687./PAR*1e2;
            rSIF760 = SIF760./PAR*1e2;
            
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
            
            
            SZA_mean = zeros(cycle_n, 1);           
            DOY_mean = zeros(cycle_n, 1);
            
            
            allsample = 1:counts;
            filters_random = randsample(allsample(SIF687>0.2 & SIF760>0.2 & RAA<180 & VZA>0),obs_numbers);
            %filters_noon = Cycle_Num == 7;% & rSIF687>0.2 & rSIF760>0.2;
            [Kparms_760, R2, RMSE] = kernelParameterRetrieval(SZA(filters_random), VZA(filters_random), RAA(filters_random), rSIF760(filters_random),'', 0);
            [Kparms_687, R2, RMSE] = kernelParameterRetrieval(SZA(filters_random), VZA(filters_random), RAA(filters_random), rSIF687(filters_random), '', 0);
            [Kparms_NIR, R2, RMSE] = kernelParameterRetrieval(SZA(filters_random), VZA(filters_random), RAA(filters_random), Refl_NIR(filters_random),'', 0);
            [Kparms_Red, R2, RMSE] = kernelParameterRetrieval(SZA(filters_random), VZA(filters_random), RAA(filters_random), Refl_Red(filters_random), '', 0);
            
            
            for cycle_i = 1:cycle_n
                
                filters = Cycle_Num == cycle_i & SIF687>0.2 & SIF760>0.2;
                filters_nadir = Cycle_Num == cycle_i & VZA == 0 & SIF687>0.2 & SIF760>0.2;
                count_n = sum(filters);
                
                
                
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
                
                %         nSIF760_true(filters) = mean(rSIF760(filters_nadir));
                %         nSIF760_NIR(filters) = mean(Refl_NIR(filters_nadir))./Refl_NIR(filters).*rSIF760(filters);
                %         nSIF760_NIRv(filters) = mean(NIRv(filters_nadir))./NIRv(filters).*rSIF760(filters);
                %         %nSIF760_EVI2(filters) = mean(EVI2(filters_nadir))./EVI2(filters).*rSIF760(filters);
                %
                %         nSIF687_true(filters) = mean(rSIF687(filters_nadir));
                %         nSIF687_Red(filters) = mean(Refl_Red(filters_nadir))./Refl_Red(filters).*rSIF687(filters);
                %         nSIF687_Redv(filters) = mean(Redv(filters_nadir))./Redv(filters).*rSIF687(filters);
                %nSIF687_EVI2(filters) = mean(EVI2(filters_nadir))./EVI2(filters).*rSIF687(filters);
                
                %% normalized by kernel-driven models
                %         switch name_i
                %             case 1
                %                 num_noon = 7;
                %             case 2
                %                 num_noon = 11;
                %             case 3
                %                 num_noon = 7;
                %         end
                c_factor = CalculateCorrectionFactor(SZA(filters), VZA(filters), RAA(filters), Kparms_760);
                nSIF760_Kernel(filters) = c_factor.*rSIF760(filters);
                
                c_factor = CalculateCorrectionFactor(SZA(filters), VZA(filters), RAA(filters), Kparms_687);
                nSIF687_Kernel(filters) = c_factor.*rSIF687(filters);
                
                
                %% accuracy statistics
                % NIR
                
            end
            
            %% calculate rMAE
            rMAE760(rand_n, 1) = nanmean(abs(nSIF760_true- rSIF760)./nSIF760_true);
            rMAE760(rand_n, 2) = nanmean(abs(nSIF760_true- nSIF760_NIR)./nSIF760_true);
            rMAE760(rand_n, 3) = nanmean(abs(nSIF760_true- nSIF760_NIRv)./nSIF760_true);
            rMAE760(rand_n, 4) = nanmean(abs(nSIF760_true- nSIF760_Kernel)./nSIF760_true);
            
            rMAE687(rand_n, 1) = nanmean(abs(nSIF687_true- rSIF687)./nSIF687_true);
            rMAE687(rand_n, 2) = nanmean(abs(nSIF687_true- nSIF687_Red)./nSIF687_true);
            rMAE687(rand_n, 3) = nanmean(abs(nSIF687_true- nSIF687_Redv)./nSIF687_true);
            rMAE687(rand_n, 4) = nanmean(abs(nSIF687_true- nSIF687_Kernel)./nSIF687_true);
        end
        %% save
        save(['results/random/RAA_0_90' num2str(obs_numbers) '_random_' nameStrings{name_i} '.mat']);
    end
end