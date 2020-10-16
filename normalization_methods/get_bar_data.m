nameStrings = {'Chickpea', 'Grass', 'Rice'};

chickpea_760 = zeros(13*4, 4);
grass_760 = zeros(13*4, 4);
rice_760 = zeros(13*4, 4);

chickpea_687 = zeros(13*4, 4);
grass_687 = zeros(13*4, 4);
rice_687 = zeros(13, 4);

for obs_numbers = 4:16

load(['results/random/RAA_0_90' num2str(obs_numbers) '_random_' nameStrings{1} '.mat']);
rMAE687(rMAE687>0.3) = nan;
rMAE760(rMAE760>0.3) = nan;
index = [(obs_numbers-3) (obs_numbers-3 + 13) (obs_numbers-3 + 13*2)  (obs_numbers-3 + 13*3)];

chickpea_687(index, 1) = obs_numbers;
chickpea_687(index, 2) = 1:4;
chickpea_687(index, 3) = nanmean(rMAE687, 1);
chickpea_687(index, 4) = nanstd(rMAE687, 1);

chickpea_760(index, 1) = obs_numbers;
chickpea_760(index, 2) = 1:4;
chickpea_760(index, 3) = nanmean(rMAE760, 1);
chickpea_760(index, 4) = nanstd(rMAE760, 1);


load(['results/random/RAA_0_90' num2str(obs_numbers) '_random_' nameStrings{2} '.mat']);
rMAE687(rMAE687>0.3) = nan;
rMAE760(rMAE760>0.3) = nan;

grass_687(index, 1) = obs_numbers;
grass_687(index, 2) = 1:4;
grass_687(index, 3) = nanmean(rMAE687, 1);
grass_687(index, 4) = nanstd(rMAE687, 1);

grass_760(index, 1) = obs_numbers;
grass_760(index, 2) = 1:4;
grass_760(index, 3) = nanmean(rMAE760, 1);
grass_760(index, 4) = nanstd(rMAE760, 1);


load(['results/random/RAA_0_90' num2str(obs_numbers) '_random_' nameStrings{3} '.mat']);
rMAE687(rMAE687>0.5) = nan;
rMAE760(rMAE760>0.5) = nan;

rice_687(index, 1) = obs_numbers;
rice_687(index, 2) = 1:4;
rice_687(index, 3) = nanmean(rMAE687, 1);
rice_687(index, 4) = nanstd(rMAE687, 1);

rice_760(index, 1) = obs_numbers;
rice_760(index, 2) = 1:4;
rice_760(index, 3) = nanmean(rMAE760, 1);
rice_760(index, 4) = nanstd(rMAE760, 1);

end

save backward_all_mean_std.mat chickpea_760 grass_760 rice_760 chickpea_687  grass_687 rice_687
