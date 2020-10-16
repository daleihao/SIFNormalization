clc;
clear all;
close all;
%% import data
names = {'corn', 'wheat'};
for tt = 1:2
    name = names{tt};
    load(['../rawdata/'  name '_2018_data.mat']);
    
    %% preprocessdata
    nir_data = double(nir_data);
    red_data = double(red_data);
    sifa_data = double(sifa_data);
    sifb_data = double(sifb_data);
    gpp_data = double(gpp_data);
    sza_data = double(sza_data);
    vza_data = double(vza_data);
    saa_data = double(saa_data);
    vaa_data = double(vaa_data);
    
    NIRs = nir_data(:,3:end);
    Reds = red_data(:,3:end);
    FarSIFs = sifa_data(:,3:end);
    RedSIFs = sifb_data(:,3:end);
    
    SZAs = sza_data(:,3:end);
    VZAs = vza_data(:,3:end);
    VAAs = vaa_data(:,3:end);
    SAAs = saa_data(:,3:end);
    RAAs = abs(VAAs-SAAs);
    RAAs(RAAs>180) = 360 - RAAs(RAAs>180);
    
    lai = gpp_data(:,8);
    LAI = repmat(lai,[1 25]);
    
    par = gpp_data(:,6);
    PARs = repmat(par, [1 25]);
    
    apar = gpp_data(:,7);
    APARs = repmat(apar, [1 25]);
    
    gpp =  gpp_data(:,5);
    GPPs =  repmat(gpp, [1 25]);
    % GPPs = APARs;
    
    doy = gpp_data(:,1);
    Doys = repmat(doy, [ 1 25]);
    
    hour = gpp_data(:,2);
    Hours = repmat(hour, [ 1 25]);
    
    ci = gpp_data(:,3);
    CIs = repmat(ci, [ 1 25]);
    
    %% calculate vegetation indices
    NDVIs = (NIRs - Reds)./(NIRs + Reds);
    NIRvs = NDVIs.*NIRs;
    Redvs = NDVIs.*NDVIs.*Reds;
    
    FPARs = APARs./PARs;
    
    LUEs = GPPs./APARs;
    
    i0s = 1-exp(-0.5.*LAI./cos(SZAs.*pi/180));
    Farw = 0.9;
    Redw = 0.1;
    
    
    doy_unique = unique(Doys(:));
    
    nadirNIRs = nan(size(RAAs));
    nadirReds = nan(size(RAAs));
    fitNIRs = nan(size(RAAs));
    fitReds = nan(size(RAAs));
    fitFarSIFs = nan(size(RAAs));
    fitRedSIFs = nan(size(RAAs));
    nadirFarSIF_ks = nan(size(RAAs));
    nadirRedSIF_ks = nan(size(RAAs));
    c_factor_far = nan(size(RAAs));
    c_factor_red = nan(size(RAAs));
    
    Rows = 1:size(RAAs, 1);
    Rows = repmat(Rows', [1 25]);
    
    
    
    for i = 1:size(doy_unique, 1)
        
        %filters = Rows == i;
        filters = Doys == doy_unique(i) & NIRs>0 & NIRs<0.9 & CIs>0.4   & SZAs<=65;
        
        %% fitting reflectance
        [Kparms_NIR, R2, RMSE, Relf_NIR_est_tmp] = kernelParameterRetrieval_refl(SZAs(filters), VZAs(filters), RAAs(filters), NIRs(filters),'', 0);
        % R2_NIR(cycle_i) = R2;
        % RMSE_NIR(cycle_i) = RMSE;
        filters = Doys == doy_unique(i) & Reds>0 & Reds<1 & CIs>0.4  & SZAs<=65;
        
        [Kparms_Red, R2, RMSE, Relf_Red_est_tmp] = kernelParameterRetrieval_refl(SZAs(filters), VZAs(filters), RAAs(filters), Reds(filters), '', 0);
        % R2_Red(cycle_i) = R2;
        % RMSE_Red(cycle_i) = RMSE;
        
        filters = Doys == doy_unique(i) & NIRs>0 & NIRs<1  & Reds>0 & Reds<1 & CIs>0.4;
        count_n = sum(filters(:));
        nadirNIRs(filters)  = CalculateRefl(SZAs(filters), zeros(count_n, 1), zeros(count_n, 1), Kparms_NIR);
        nadirReds(filters)  = CalculateRefl(SZAs(filters), zeros(count_n, 1), zeros(count_n, 1), Kparms_Red);
        
        fitNIRs(filters)  = CalculateRefl(SZAs(filters), VZAs(filters), RAAs(filters), Kparms_NIR);
        fitReds(filters)  = CalculateRefl(SZAs(filters), VZAs(filters), RAAs(filters), Kparms_Red);
        
        %% fitting SIF
        filters = Doys == doy_unique(i) & FarSIFs>0 & FarSIFs<4 & PARs>0  & CIs>0.4  & SZAs<=65;
        
        [Kparms, R2, RMSE, SIF760_est_tmp] = kernelParameterRetrieval_refl(SZAs(filters), VZAs(filters), RAAs(filters), FarSIFs(filters)./PARs(filters), '', 0);
        filters = Doys == doy_unique(i) & FarSIFs>0 & FarSIFs<4 & PARs>0  & CIs>0.4;
        c_factor_far(filters) = CalculateCorrectionFactor(SZAs(filters), VZAs(filters), RAAs(filters), Kparms);
        nadirFarSIF_ks(filters) = c_factor_far(filters).*FarSIFs(filters);
        fitFarSIFs(filters) = CalculateRefl(SZAs(filters), VZAs(filters), RAAs(filters), Kparms).*PARs(filters);
        %R2_760(cycle_i) = R2;
        %RMSE_760(cycle_i) = RMSE;
        filters = Doys == doy_unique(i) & RedSIFs>0 & RedSIFs<4  & PARs>0  & CIs>0.4 & SZAs<=65;
        
        [Kparms, R2, RMSE, SIF687_est_tmp] = kernelParameterRetrieval_refl(SZAs(filters), VZAs(filters), RAAs(filters), RedSIFs(filters)./PARs(filters), '', 0);
        % CalculateCorrectionFactor(SZAs(filters), VZAs(filters), RAAs(filters), Kparms)
        filters = Doys == doy_unique(i) & RedSIFs>0 & RedSIFs<4  & PARs>0  & CIs>0.4 ;
        c_factor_red(filters) = CalculateCorrectionFactor(SZAs(filters), VZAs(filters), RAAs(filters), Kparms);
        nadirRedSIF_ks(filters) = c_factor_red(filters).*RedSIFs(filters);
        fitRedSIFs(filters) = CalculateRefl(SZAs(filters), VZAs(filters), RAAs(filters), Kparms).*PARs(filters);
        
        %R2_687(cycle_i) = R2;
        %RMSE_687(cycle_i) = RMSE;
        %kRed0(filters)  = CalculateRefl(SZA(filters), zeros(count_n, 1), zeros(count_n, 1), Kparms_Red);
    end
    
    
    nadirNDVIs = (nadirNIRs-nadirReds)./(nadirNIRs-nadirReds);
    nadirNIRvs = nadirNIRs.*nadirNDVIs;
    nadirRedvs = nadirReds.*nadirNDVIs.*nadirNDVIs;
    
    fitNDVIs = (fitNIRs-fitReds)./(fitNIRs-fitReds);
    fitNIRvs = fitNIRs.*fitNDVIs;
    fitRedvs = fitReds.*fitNDVIs.*fitNDVIs;
    
    
    
    FarFescs = NIRvs./FPARs./Farw;
    RedFescs = Redvs./FPARs./Redw;
    
    nadirFarFescs = nadirNIRvs./FPARs./Farw;
    nadirRedFescs = nadirRedvs./FPARs./Redw;
    nadirFarFescs_i0 = nadirNIRvs./i0s./Farw;
    nadirRedFescs_i0 = nadirRedvs./i0s./Redw;
    
    %% far-red
    nadirFarSIF_rs = FarSIFs./NIRvs.*nadirNIRvs;
    totalFarSIF_fpars = FarSIFs./FarFescs;
    totalFarSIF_i0s = FarSIFs./NIRvs.*i0s.*Farw;
    %% red
    nadirRedSIF_rs = RedSIFs./Redvs.*nadirRedvs;
    totalRedSIF_fpars = RedSIFs./RedFescs;
    totalRedSIF_i0s = RedSIFs./Redvs.*i0s.*Redw;
    
    
    FarSIFyields =  nanmean(nadirRedSIF_rs./APARs,2);
    FarPhiFs = nanmean(totalRedSIF_fpars./APARs,2);
    
    RedSIFyields =  RedSIFs./APARs;
    RedPhiFs = RedSIFs./APARs./RedFescs;
    
    
    
    %     figure;
    %     subplot(241)
    %     plot(GPPs, 'r')
    %     subplot(242)
    %       plot(FarSIFs(:,25), 'r')
    %      subplot(243)
    %     plot(nadirFarSIF_ks(:,25), 'r')
    %     subplot(244)
    %     plot(nadirFarSIF_rs(:,25), 'r')
    %     subplot(245)
    %     plot(LUEs, 'r')
    %      subplot(246)
    %     plot(nadirNIRvs./FPARs, 'r')
    %      subplot(247)
    %     plot(FPARs, 'r')
    %%
    Cols = 1:size(RAAs, 2);
    Cols = repmat(Cols, [size(RAAs, 1) 1]);
    R2s_far = zeros(25, 7);
    R2s_red = zeros(25, 7);
    
    %%
    
    for col_num = 1:25
        filters =  GPPs>0 & ...
            FarSIFs>0 & nadirFarSIF_rs>0  & nadirFarSIF_ks>0 & totalFarSIF_fpars>0 & totalFarSIF_i0s>0 & ...
            FarSIFs<30 & nadirFarSIF_rs<30 & nadirFarSIF_ks<30 & totalFarSIF_fpars<30 & totalFarSIF_i0s<30 & ...
            Cols == col_num & FarFescs<3 & LUEs>0.01 & ...
            c_factor_far>0 & c_factor_far<4  & APARs>0  & CIs>0.4;
        
        GPP_filters = GPPs(filters);
        FarSIF_filters = FarSIFs(filters);
        nadirFarSIF_r_filters = nadirFarSIF_rs(filters);
        nadirFarSIF_k_filters = nadirFarSIF_ks(filters);
        totalFarSIF_fpar_filters = totalFarSIF_fpars(filters);
        totalFarSIF_i0_filters = totalFarSIF_i0s(filters);
        APAR_filters = APARs(filters);
        APAR_Fesc_filters = FarFescs(filters);
        APAR_Fesc2_filters = nadirFarFescs(filters);
        LUE_filters = LUEs(filters);
        
        startPoints = [0 0];
        
        mdl1 = HyperbolicFit(FarSIF_filters, GPP_filters,startPoints);
        mdl2 = HyperbolicFit(nadirFarSIF_r_filters, GPP_filters,startPoints);
        mdl3 = HyperbolicFit(nadirFarSIF_k_filters, GPP_filters,startPoints);
        mdl4 = HyperbolicFit(totalFarSIF_fpar_filters, GPP_filters,startPoints);
        mdl5 = HyperbolicFit(totalFarSIF_i0_filters, GPP_filters,startPoints);
        mdl6 = HyperbolicFit(APAR_filters, GPP_filters,startPoints);
        mdl7 = HyperbolicFit(APAR_Fesc_filters, LUE_filters,startPoints);
        mdl8 = HyperbolicFit(APAR_Fesc2_filters, LUE_filters,startPoints);
        
        R2s_far(col_num, 1) = mdl1.Rsquared.Ordinary;
        R2s_far(col_num, 2) = mdl2.Rsquared.Ordinary;
        R2s_far(col_num, 3) = mdl3.Rsquared.Ordinary;
        R2s_far(col_num, 4) = mdl4.Rsquared.Ordinary;
        R2s_far(col_num, 5) = mdl5.Rsquared.Ordinary;
        R2s_far(col_num, 6) = mdl6.Rsquared.Ordinary;
        R2s_far(col_num, 7) = mdl7.Rsquared.Ordinary;
        R2s_far(col_num, 8) = mdl8.Rsquared.Ordinary;
        
        %     figure;
        %     hold on
        %     scatter(FarSIF_filters, GPP_filters, 'r')
        %     scatter(nadirFarSIF_r_filters, GPP_filters, 'b')
        %     scatter(nadirFarSIF_k_filters, GPP_filters,'k')
        
        %         startPoints = [0.6 0.6];
        %         [fitPsh1, staths1] = HyperbolicFit(FarSIFs, GPPs, startPoints);
        %         [fitPsh2, staths2] = HyperbolicFit(nadirFarSIF_rs, GPPs, startPoints);
        %         [fitPsh3, staths3] = HyperbolicFit(nadirFarSIF_ks, GPPs, startPoints);
        %         [fitPsh4, staths4] = HyperbolicFit(totalFarSIF_fpars, GPPs, startPoints);
        %         [fitPsh5, staths5] = HyperbolicFit(totalFarSIF_fpars, GPPs, startPoints);
        %
        
        
        filters =  GPPs>0 & ...
            RedSIFs>0 & nadirRedSIF_rs>0  & nadirRedSIF_ks>0 & totalRedSIF_fpars>0 & totalRedSIF_i0s>0 & ...
            RedSIFs<30 & nadirRedSIF_rs<30 & nadirRedSIF_ks<30 & totalRedSIF_fpars<30 & totalRedSIF_i0s<30 & ...
            Cols == col_num  & RedFescs<3 & LUEs>0.01 & ...
            c_factor_red>0 & c_factor_red<4 & APARs>0  & CIs>0.4;
        % disp(sum(filters(:)))
        GPP_filters = GPPs(filters);
        RedSIF_filters = RedSIFs(filters);
        nadirRedSIF_r_filters = nadirRedSIF_rs(filters);
        nadirRedSIF_k_filters = nadirRedSIF_ks(filters);
        totalRedSIF_fpar_filters = totalRedSIF_fpars(filters);
        totalRedSIF_i0_filters = totalRedSIF_i0s(filters);
        APAR_filters = APARs(filters);
        APAR_Fesc_filters = RedFescs(filters);
        APAR_Fesc2_filters = nadirRedFescs(filters);
        LUE_filters = LUEs(filters);
        
        startPoints = [0 0];
        
        mdl1 = HyperbolicFit(RedSIF_filters, GPP_filters,startPoints);
        mdl2 = HyperbolicFit(nadirRedSIF_r_filters, GPP_filters,startPoints);
        mdl3 = HyperbolicFit(nadirRedSIF_k_filters, GPP_filters,startPoints);
        mdl4 = HyperbolicFit(totalRedSIF_fpar_filters, GPP_filters,startPoints);
        mdl5 = HyperbolicFit(totalRedSIF_i0_filters, GPP_filters,startPoints);
        mdl6 = HyperbolicFit(APAR_filters, GPP_filters,startPoints);
        mdl7 = HyperbolicFit(APAR_Fesc_filters, LUE_filters,startPoints);
        mdl8 = HyperbolicFit(APAR_Fesc2_filters, LUE_filters,startPoints);
        
        
        
        
        R2s_red(col_num, 1) = mdl1.Rsquared.Ordinary;
        R2s_red(col_num, 2) = mdl2.Rsquared.Ordinary;
        R2s_red(col_num, 3) = mdl3.Rsquared.Ordinary;
        R2s_red(col_num, 4) = mdl4.Rsquared.Ordinary;
        R2s_red(col_num, 5) = mdl5.Rsquared.Ordinary;
        R2s_red(col_num, 6) = mdl6.Rsquared.Ordinary;
        R2s_red(col_num, 7) = mdl7.Rsquared.Ordinary;
        R2s_red(col_num, 8) = mdl8.Rsquared.Ordinary;
        %     figure;
        %     hold on
        %     scatter(RedSIF_filters, GPP_filters, 'r')
        %     scatter(nadirRedSIF_r_filters, GPP_filters, 'b')
        %     scatter(nadirRedSIF_k_filters, GPP_filters,'k')
        %         startPoints = [0.1 1];
        %
        %         [fitPsh1, staths1] = HyperbolicFit(RedSIF_filters, GPP_filters, startPoints);
        %
        %         startPoints = [10 0];
        %         [fitPsh2, staths2] = HyperbolicFit(nadirRedSIF_r_filters, GPP_filters, startPoints);
        %         startPoints = [0.6 2];
        %
        %         [fitPsh3, staths3] = HyperbolicFit(nadirRedSIF_k_filters, GPP_filters, startPoints);
        %         [fitPsh4, staths4] = HyperbolicFit(totalRedSIF_fpar_filters, GPP_filters, startPoints);
        %         [fitPsh5, staths5] = HyperbolicFit(totalRedSIF_fpar_filters, GPP_filters, startPoints);
    end
    
    save(['results/hyper_R2s_' name '.mat'], 'R2s_far', 'R2s_red');
    
end

