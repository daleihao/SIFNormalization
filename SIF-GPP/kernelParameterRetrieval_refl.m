function [Kparms, R2, RMSE, Refl_e] = kernelParameterRetrieval_refl(SZA, VZA, RAA, Refl, filename, isCal)

%% deg to rad

SZA_0 = SZA;
VZA_0 = VZA;
RAA_0 = RAA;

SZA = abs(SZA)*pi/180;
VZA = abs(VZA)*pi/180;
RAA = abs(RAA)*pi/180;

count = size(SZA, 1);
Kiso = ones(count, 1);
Kvol = vkeroThick(SZA, VZA, RAA);
Kgeo = LiTransit(SZA, VZA, RAA);
Kparms= [Kiso Kvol Kgeo]\Refl;

Refl_e = [Kiso Kvol Kgeo]*Kparms;

if size(Refl_e, 1)>1 | size(Refl_e, 2)>1 
    R = corrcoef(Refl_e, Refl);
    R = R(1,2);
    R2 = R^2;
else
    R2 = nan;
end
RMSE = sqrt(mean((Refl_e-Refl).^2));



%% plot PP and SPP
if(isCal)
    VZA_0(RAA_0==0) = -VZA_0(RAA_0==0);
    VZA_0(RAA_0==90) = -VZA_0(RAA_0==90);
    figure;
    subplot(1,2,1)
    hold on
    plot(VZA_0(RAA_0==0 | RAA_0==180), Refl(RAA_0==0 | RAA_0==180), 'ro')
    plot(VZA_0(RAA_0==0 | RAA_0==180), Refl_e(RAA_0==0 | RAA_0==180), 'b*')
    title(['PP' num2str(mean(SZA_0), '%.2f')])
    subplot(1,2,2)
    hold on
    plot(VZA_0(RAA_0==90 | RAA_0==270), Refl(RAA_0==90 | RAA_0==270), 'ro')
    plot(VZA_0(RAA_0==90 | RAA_0==270), Refl_e(RAA_0==90 | RAA_0==270), 'b*')
    title(['CPP: R2 = ' num2str(R2, '%.2f') ' RMSE = ' num2str(RMSE, '%.2f')])
    print(gcf, '-dtiff', '-r600',['results/SIF_fit/' filename '.tif']);
    close all;
end
end