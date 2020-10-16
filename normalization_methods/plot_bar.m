nameStrings = {'Chickpea', 'Grass', 'Rice'};
chickpea_mean_760 = zeros(13, 4);
grass_mean_760 = zeros(13, 4);
rice_mean_760 = zeros(13, 4);
chickpea_std_760 = zeros(13, 4);
grass_std_760 = zeros(13, 4);
rice_std_760 = zeros(13, 4);

chickpea_mean_687 = zeros(13, 4);
grass_mean_687 = zeros(13, 4);
rice_mean_687 = zeros(13, 4);
chickpea_std_687 = zeros(13, 4);
grass_std_687 = zeros(13, 4);
rice_std_687 = zeros(13, 4);
for obs_numbers = 4:16

load(['results/random/' num2str(obs_numbers) '_random_' nameStrings{1} '.mat']);
rMAE687(rMAE687>0.3) = nan;
rMAE760(rMAE760>0.3) = nan;
chickpea_mean_687(obs_numbers-3, :) = nanmean(rMAE687, 1);
chickpea_mean_687(obs_numbers-3, :) = nanmean(rMAE687, 1);
chickpea_std_687(obs_numbers-3, :) = nanstd(rMAE687, 1);
chickpea_std_687(obs_numbers-3, :) = nanstd(rMAE687, 1);
chickpea_mean_760(obs_numbers-3, :) = nanmean(rMAE760, 1);
chickpea_mean_760(obs_numbers-3, :) = nanmean(rMAE760, 1);
chickpea_std_760(obs_numbers-3, :) = nanstd(rMAE760, 1);
chickpea_std_760(obs_numbers-3, :) = nanstd(rMAE760, 1);

load(['results/random/' num2str(obs_numbers) '_random_' nameStrings{2} '.mat']);
rMAE687(rMAE687>0.3) = nan;
rMAE760(rMAE760>0.3) = nan;
grass_mean_687(obs_numbers-3, :) = nanmean(rMAE687, 1);
grass_mean_687(obs_numbers-3, :) = nanmean(rMAE687, 1);
grass_std_687(obs_numbers-3, :) = nanstd(rMAE687, 1);
grass_std_687(obs_numbers-3, :) = nanstd(rMAE687, 1);
grass_mean_760(obs_numbers-3, :) = nanmean(rMAE760, 1);
grass_mean_760(obs_numbers-3, :) = nanmean(rMAE760, 1);
grass_std_760(obs_numbers-3, :) = nanstd(rMAE760, 1);
grass_std_760(obs_numbers-3, :) = nanstd(rMAE760, 1);

load(['results/random/' num2str(obs_numbers) '_random_' nameStrings{3} '.mat']);
rMAE687(rMAE687>0.5) = nan;
rMAE760(rMAE760>0.5) = nan;
rice_mean_687(obs_numbers-3, :) = nanmean(rMAE687, 1);
rice_mean_687(obs_numbers-3, :) = nanmean(rMAE687, 1);
rice_std_687(obs_numbers-3, :) = nanstd(rMAE687, 1);
rice_std_687(obs_numbers-3, :) = nanstd(rMAE687, 1);
rice_mean_760(obs_numbers-3, :) = nanmean(rMAE760, 1);
rice_mean_760(obs_numbers-3, :) = nanmean(rMAE760, 1);
rice_std_760(obs_numbers-3, :) = nanstd(rMAE760, 1);
rice_std_760(obs_numbers-3, :) = nanstd(rMAE760, 1);

end

