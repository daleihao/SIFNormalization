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
    
    FarFescs = NIRs./FPARs;
    RedFescs = Reds./FPARs;
    
    nadirFarFescs = nadirNIRvs./FPARs;
    nadirRedFescs = nadirRedvs./FPARs;
    %% far-red
    nadirFarSIF_rs = FarSIFs./NIRs.*nadirNIRs;
    totalFarSIF_fpars = FarSIFs./FarFescs;
    totalFarSIF_i0s = FarSIFs./NIRs.*i0s.*Farw;
    %% red
    nadirRedSIF_rs = RedSIFs./Reds.*nadirReds;
    totalRedSIF_fpars = RedSIFs./RedFescs;
    totalRedSIF_i0s = RedSIFs./Reds.*i0s.*Redw;
    
    
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
        
        [fitPs1, stats1] = LinearFit(FarSIF_filters, GPP_filters);
        [fitPs2, stats2] = LinearFit(nadirFarSIF_r_filters, GPP_filters);
        [fitPs3, stats3] = LinearFit(nadirFarSIF_k_filters, GPP_filters);
        [fitPs4, stats4] = LinearFit(totalFarSIF_fpar_filters, GPP_filters);
        [fitPs5, stats5] = LinearFit(totalFarSIF_i0_filters, GPP_filters);
        [fitPs6, stats6] = LinearFit(APAR_filters, GPP_filters);
        [fitPs7, stats7] = LinearFit(APAR_Fesc_filters, LUE_filters);
        [fitPs8, stats8] = LinearFit(APAR_Fesc2_filters, LUE_filters);
        
        R2s_far(col_num, 1) = stats1.rsquare;
        R2s_far(col_num, 2) = stats2.rsquare;
        R2s_far(col_num, 3) = stats3.rsquare;
        R2s_far(col_num, 4) = stats4.rsquare;
        R2s_far(col_num, 5) = stats5.rsquare;
        R2s_far(col_num, 6) = stats6.rsquare;
        R2s_far(col_num, 7) = stats7.rsquare;
        R2s_far(col_num, 8) = stats8.rsquare;
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
        
        
        [fitPs1, stats1] = LinearFit(RedSIF_filters, GPP_filters);
        [fitPs2, stats2] = LinearFit(nadirRedSIF_r_filters, GPP_filters);
        [fitPs3, stats3] = LinearFit(nadirRedSIF_k_filters, GPP_filters);
        [fitPs4, stats4] = LinearFit(totalRedSIF_fpar_filters, GPP_filters);
        [fitPs5, stats5] = LinearFit(totalRedSIF_i0_filters, GPP_filters);
        [fitPs6, stats6] = LinearFit(APAR_filters, GPP_filters);
        [fitPs7, stats7] = LinearFit(APAR_Fesc_filters, LUE_filters);
        [fitPs8, stats8] = LinearFit(APAR_Fesc2_filters, LUE_filters);
        
        R2s_red(col_num, 1) = stats1.rsquare;
        R2s_red(col_num, 2) = stats2.rsquare;
        R2s_red(col_num, 3) = stats3.rsquare;
        R2s_red(col_num, 4) = stats4.rsquare;
        R2s_red(col_num, 5) = stats5.rsquare;
        R2s_red(col_num, 6) = stats6.rsquare;
        R2s_red(col_num, 7) = stats7.rsquare;
        R2s_red(col_num, 8) = stats8.rsquare;
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
    
    
    Rows = 1:size(Rows, 1);
    Rows = repmat(Rows', [1 size(RAAs, 2)]);
    
    CVs_far = zeros(size(Rows, 1), 5);
    CVs_red = zeros(size(Rows, 1), 5);
    %%
    
    for row_num = 1:size(Rows, 1)
        filters =  GPPs>0 & ...
            FarSIFs>0 & nadirFarSIF_rs>0  & nadirFarSIF_ks>0 & totalFarSIF_fpars>0 & totalFarSIF_i0s>0 & ...
            FarSIFs<30 & nadirFarSIF_rs<30 & nadirFarSIF_ks<30 & totalFarSIF_fpars<30 & totalFarSIF_i0s<30 & ...
            Rows == row_num & FarFescs<3 & LUEs>0.01 & ...
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
        
        
        CVs_far(row_num, 1) = nanstd(FarSIF_filters)/nanmean(FarSIF_filters);
        CVs_far(row_num, 2) = nanstd(nadirFarSIF_r_filters)/nanmean(nadirFarSIF_r_filters);
        CVs_far(row_num, 3) = nanstd(nadirFarSIF_k_filters)/nanmean(nadirFarSIF_k_filters);
        CVs_far(row_num, 4) = nanstd(totalFarSIF_fpar_filters)/nanmean(totalFarSIF_fpar_filters);
        CVs_far(row_num, 5) = nanstd(totalFarSIF_i0_filters)/nanmean(totalFarSIF_i0_filters);
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
        
        
        filters = GPPs>0 & ...
            RedSIFs>0 & nadirRedSIF_rs>0  & nadirRedSIF_ks>0 & totalRedSIF_fpars>0 & totalRedSIF_i0s>0 & ...
            RedSIFs<30 & nadirRedSIF_rs<30 & nadirRedSIF_ks<30 & totalRedSIF_fpars<30 & totalRedSIF_i0s<30 & ...
           Rows == row_num & RedFescs<3 & LUEs>0.01 & ...
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
        
        
        CVs_red(row_num, 1) = nanstd(RedSIF_filters)/nanmean(RedSIF_filters);
        CVs_red(row_num, 2) = nanstd(nadirRedSIF_r_filters)/nanmean(nadirRedSIF_r_filters);
        CVs_red(row_num, 3) = nanstd(nadirRedSIF_k_filters)/nanmean(nadirRedSIF_k_filters);
        CVs_red(row_num, 4) = nanstd(totalRedSIF_fpar_filters)/nanmean(totalRedSIF_fpar_filters);
        CVs_red(row_num, 5) = nanstd(totalRedSIF_i0_filters)/nanmean(totalRedSIF_i0_filters);
        
    end
    %%
    figure;
    set(gcf,'unit','normalized','position',[0.1,0.1,0.6,0.6]);
    set(gca, 'Position', [0 0 1 1])
    
    subplot('Position', [ 0.02 0.1 0.45 0.8])
    
    boxplot(R2s_far(:,1:5), 'widths', 0.5)
    ylabel('R^2')
    ylim([0 1])
    set(gca, 'xticklabel', {'raw data', 'nadir_Rfl', 'nadir_KD', 'total_FPAR', 'total_i0'},...%,'APAR','fesc', 'fesc2'},...
        'linewidth', 1)
    title('Far-red')
    
    subplot('Position', [ 0.55 0.1 0.4 0.8])
    boxplot(R2s_red(:,1:5))
    %ylabel('R^2')
    ylim([0 1])
    set(gca, 'xticklabel', {'raw data', 'nadir_Rfl', 'nadir_KD', 'total_FPAR', 'total_i0'},...%,'APAR','fesc', 'fesc2'},...
        'linewidth', 1)
    title('Red')
    suptitle(name)
 
    
      %%
    figure;
    set(gcf,'unit','normalized','position',[0.1,0.1,0.6,0.6]);
    set(gca, 'Position', [0 0 1 1])
    
    subplot('Position', [ 0.02 0.1 0.45 0.8])
    
    boxplot(CVs_far(:,1:5), 'widths', 0.5)
    ylabel('CV')
    ylim([0 1])
    set(gca, 'xticklabel', {'raw data', 'nadir_Rfl', 'nadir_KD', 'total_FPAR', 'total_i0'},...
        'linewidth', 1)
    title('Far-red')
    
    subplot('Position', [ 0.55 0.1 0.4 0.8])
    boxplot(CVs_red(:,1:5))
    %ylabel('R^2')
    ylim([0 1])
    set(gca, 'xticklabel', {'raw data', 'nadir_Rfl', 'nadir_KD', 'total_FPAR', 'total_i0'},...
        'linewidth', 1)
    title('Red')
    suptitle(name)
    
    save([ 'results/CVs_' name '_nir.mat'], 'CVs_far', 'CVs_red');
    save(['results/R2s_' name '_nir.mat'], 'R2s_far', 'R2s_red');
    %%
    %     figure;
    %     subplot(121)
    %     boxplot(R2s_far-R2s_far(:,1))
    %     ylabel('delta R^2')
    %     set(gca, 'xticklabel', {'raw data', 'nadir_Rfl', 'nadir_KD', 'total_FPAR', 'total_i0','APAR'},...
    %         'linewidth', 1)
    %
    %     subplot(122)
    %     boxplot(R2s_red-R2s_red(:,1))
    %     ylabel('delta R^2')
    %     set(gca, 'xticklabel', {'raw data', 'nadir_Rfl', 'nadir_KD', 'total_FPAR', 'total_i0','APAR'},...
    %         'linewidth', 1)
    %     suptitle(name)
    
    
    
    
    %     %% SAA mean
    %     R2s_far = zeros(1, 7);
    %     R2s_red = zeros(1, 7);
    %     col_num = 1;
    %     filters = GPPs>0 & ...
    %         FarSIFs>0 & nadirFarSIF_rs>0  & nadirFarSIF_ks>0 & totalFarSIF_fpars>0 & totalFarSIF_i0s>0 & ...
    %         FarSIFs<30 & nadirFarSIF_rs<30 & nadirFarSIF_ks<30 & totalFarSIF_fpars<30 & totalFarSIF_i0s<30 & ...
    %         FarFescs<1 & LUEs>0.01 & ...
    %         c_factor_far>0 & c_factor_far<4  & APARs>0;
    %
    %     GPP_alls = GPPs;
    %     FarSIF_alls = FarSIFs;
    %     nadirFarSIF_r_alls = nadirFarSIF_rs;
    %     nadirFarSIF_k_alls = nadirFarSIF_ks;
    %     totalFarSIF_fpar_alls = totalFarSIF_fpars;
    %     totalFarSIF_i0_alls = totalFarSIF_i0s;
    %     APAR_alls = APARs;
    %     APAR_Fesc_alls = APARs.*FarFescs;
    %
    %
    %     GPP_alls(~filters) = nan;
    %     FarSIF_alls(~filters) = nan;
    %     nadirFarSIF_r_alls(~filters) = nan;
    %     nadirFarSIF_k_alls(~filters) = nan;
    %     totalFarSIF_fpar_alls(~filters) = nan;
    %     totalFarSIF_i0_alls(~filters) = nan;
    %     APAR_alls(~filters) = nan;
    %     APAR_Fesc_alls(~filters) = nan;
    %
    %
    %     GPP_filters = nanmean(GPP_alls,2);
    %     FarSIF_filters = nanmean(FarSIF_alls,2);
    %     nadirFarSIF_r_filters = nanmean(nadirFarSIF_r_alls,2);
    %     nadirFarSIF_k_filters = nanmean(nadirFarSIF_k_alls,2);
    %     totalFarSIF_fpar_filters = nanmean(totalFarSIF_fpar_alls, 2);
    %     totalFarSIF_i0_filters = nanmean(totalFarSIF_i0_alls, 2);
    %     APAR_filters = nanmean(APAR_alls, 2);
    %     APAR_Fesc_filters = nanmean(APAR_Fesc_alls, 2);
    %
    %     [fitPs1, stats1] = LinearFit(FarSIF_filters, GPP_filters);
    %     [fitPs2, stats2] = LinearFit(nadirFarSIF_r_filters, GPP_filters);
    %     [fitPs3, stats3] = LinearFit(nadirFarSIF_k_filters, GPP_filters);
    %     [fitPs4, stats4] = LinearFit(totalFarSIF_fpar_filters, GPP_filters);
    %     [fitPs5, stats5] = LinearFit(totalFarSIF_fpar_filters, GPP_filters);
    %     [fitPs6, stats6] = LinearFit(APAR_filters, GPP_filters);
    %     [fitPs7, stats7] = LinearFit(APAR_Fesc_filters, GPP_filters);
    %
    %     R2s_far(col_num, 1) = stats1.rsquare;
    %     R2s_far(col_num, 2) = stats2.rsquare;
    %     R2s_far(col_num, 3) = stats3.rsquare;
    %     R2s_far(col_num, 4) = stats4.rsquare;
    %     R2s_far(col_num, 5) = stats5.rsquare;
    %     R2s_far(col_num, 6) = stats6.rsquare;
    %     R2s_far(col_num, 7) = stats7.rsquare;
    %     %     figure;
    %     %     hold on
    %     %     scatter(FarSIF_filters, GPP_filters, 'r')
    %     %     scatter(nadirFarSIF_r_filters, GPP_filters, 'b')
    %     %     scatter(nadirFarSIF_k_filters, GPP_filters,'k')
    %
    %     startPoints = [0.6 0.6];
    %     [fitPsh1, staths1] = HyperbolicFit(FarSIFs, GPPs, startPoints);
    %     [fitPsh2, staths2] = HyperbolicFit(nadirFarSIF_rs, GPPs, startPoints);
    %     [fitPsh3, staths3] = HyperbolicFit(nadirFarSIF_ks, GPPs, startPoints);
    %     [fitPsh4, staths4] = HyperbolicFit(totalFarSIF_fpars, GPPs, startPoints);
    %     [fitPsh5, staths5] = HyperbolicFit(totalFarSIF_fpars, GPPs, startPoints);
    %
    %
    %     %% reds
    %     filters = GPPs>0 & ...
    %         RedSIFs>0 & nadirRedSIF_rs>0  & nadirRedSIF_ks>0 & totalRedSIF_fpars>0 & totalRedSIF_i0s>0 & ...
    %         RedSIFs<30 & nadirRedSIF_rs<30 & nadirRedSIF_ks<30 & totalRedSIF_fpars<30 & totalRedSIF_i0s<30 & ...
    %         RedFescs<1 & LUEs>0.01 & ...
    %         c_factor_red>0 & c_factor_red<4 & APARs>0;
    %
    %     GPP_alls = GPPs;
    %     RedSIF_alls = RedSIFs;
    %     nadirRedSIF_r_alls = nadirRedSIF_rs;
    %     nadirRedSIF_k_alls = nadirRedSIF_ks;
    %     totalRedSIF_fpar_alls = totalRedSIF_fpars;
    %     totalRedSIF_i0_alls = totalRedSIF_i0s;
    %     APAR_alls = APARs;
    %     APAR_Fesc_alls = APARs.*RedFescs;
    %
    %
    %     GPP_alls(~filters) = nan;
    %     RedSIF_alls(~filters) = nan;
    %     nadirRedSIF_r_alls(~filters) = nan;
    %     nadirRedSIF_k_alls(~filters) = nan;
    %     totalRedSIF_fpar_alls(~filters) = nan;
    %     totalRedSIF_i0_alls(~filters) = nan;
    %     APAR_alls(~filters) = nan;
    %     APAR_Fesc_alls(~filters) = nan;
    %
    %
    %     GPP_filters = nanmean(GPP_alls,2);
    %     RedSIF_filters = nanmean(RedSIF_alls,2);
    %     nadirRedSIF_r_filters = nanmean(nadirRedSIF_r_alls,2);
    %     nadirRedSIF_k_filters = nanmean(nadirRedSIF_k_alls,2);
    %     totalRedSIF_fpar_filters = nanmean(totalRedSIF_fpar_alls, 2);
    %     totalRedSIF_i0_filters = nanmean(totalRedSIF_i0_alls, 2);
    %     APAR_filters = nanmean(APAR_alls, 2);
    %     APAR_Fesc_filters = nanmean(APAR_Fesc_alls, 2);
    %
    %     [fitPs1, stats1] = LinearFit(RedSIF_filters, GPP_filters);
    %     [fitPs2, stats2] = LinearFit(nadirRedSIF_r_filters, GPP_filters);
    %     [fitPs3, stats3] = LinearFit(nadirRedSIF_k_filters, GPP_filters);
    %     [fitPs4, stats4] = LinearFit(totalRedSIF_fpar_filters, GPP_filters);
    %     [fitPs5, stats5] = LinearFit(totalRedSIF_fpar_filters, GPP_filters);
    %     [fitPs6, stats6] = LinearFit(APAR_filters, GPP_filters);
    %     [fitPs7, stats7] = LinearFit(APAR_Fesc_filters, GPP_filters);
    %
    %     R2s_red(col_num, 1) = stats1.rsquare;
    %     R2s_red(col_num, 2) = stats2.rsquare;
    %     R2s_red(col_num, 3) = stats3.rsquare;
    %     R2s_red(col_num, 4) = stats4.rsquare;
    %     R2s_red(col_num, 5) = stats5.rsquare;
    %     R2s_red(col_num, 6) = stats6.rsquare;
    %     R2s_red(col_num, 7) = stats7.rsquare;
    %     %     figure;
    %     %     hold on
    %     %     scatter(RedSIF_filters, GPP_filters, 'r')
    %     %     scatter(nadirRedSIF_r_filters, GPP_filters, 'b')
    %     %     scatter(nadirRedSIF_k_filters, GPP_filters,'k')
    %     startPoints = [0.1 1];
    %
    %     [fitPsh1, staths1] = HyperbolicFit(RedSIF_filters, GPP_filters, startPoints);
    %
    %     startPoints = [10 0];
    %     [fitPsh2, staths2] = HyperbolicFit(nadirRedSIF_r_filters, GPP_filters, startPoints);
    %     startPoints = [0.6 2];
    %
    %     [fitPsh3, staths3] = HyperbolicFit(nadirRedSIF_k_filters, GPP_filters, startPoints);
    %     [fitPsh4, staths4] = HyperbolicFit(totalRedSIF_fpar_filters, GPP_filters, startPoints);
    %     [fitPsh5, staths5] = HyperbolicFit(totalRedSIF_fpar_filters, GPP_filters, startPoints);
    %
    %
    %     %%
    %     figure;
    %     subplot(121)
    %     bar(R2s_far)
    %     ylabel('R^2')
    %     set(gca, 'xticklabel', {'raw data', 'nadir_Rfl', 'nadir_KD', 'total_FPAR', 'total_i0','APAR','APAR*fesc'},...
    %         'linewidth', 1)
    %
    %     subplot(122)
    %     bar(R2s_red)
    %     ylabel('R^2')
    %     set(gca, 'xticklabel', {'raw data', 'nadir_Rfl', 'nadir_KD', 'total_FPAR', 'total_i0','APAR','APAR*fesc'},...
    %         'linewidth', 1)
    %     suptitle(name)
    %
    %     %%
    %     figure;
    %     subplot(121)
    %     bar(R2s_far-R2s_far(:,1))
    %     ylabel('delta R^2')
    %     set(gca, 'xticklabel', {'raw data', 'nadir_Rfl', 'nadir_KD', 'total_FPAR', 'total_i0','APAR'},...
    %         'linewidth', 1)
    %
    %     subplot(122)
    %     bar(R2s_red-R2s_red(:,1))
    %     ylabel('delta R^2')
    %     set(gca, 'xticklabel', {'raw data', 'nadir_Rfl', 'nadir_KD', 'total_FPAR', 'total_i0','APAR'},...
    %         'linewidth', 1)
    %     suptitle(name)
    
    
end

% figure;
% hold on
%
% scatter(nadir_farsif(filters), gpp(filters),'r')
% scatter(total_farsif(filters), gpp(filters),'g')
% scatter(total_farsif2(filters), gpp(filters),'k')
%
% nadir_redsif = redsif./red.*nadir_red;
% total_redsif = redsif./red.*fpar;
% total_redsif2 = redsif./red.*i0.*0.1;
% filters = redsif>0 & nadir_redsif>0 & gpp>0 & total_redsif>0 & total_redsif2>0 & nadir_redsif<100 & total_redsif<100 & total_redsif2<100  & 239<VAA & VAA<241;
% disp(corrcoef(redsif(filters), gpp(filters)))
% disp(corrcoef(nadir_redsif(filters), gpp(filters)))
% disp(corrcoef(total_redsif(filters), gpp(filters)))
% disp(corrcoef(total_redsif2(filters), gpp(filters)))
%
% figure;
% hold on
%
% scatter(redsif(filters), gpp(filters),'r')
% scatter(nadir_redsif(filters), gpp(filters),'g')
% scatter(total_redsif(filters), gpp(filters),'k')