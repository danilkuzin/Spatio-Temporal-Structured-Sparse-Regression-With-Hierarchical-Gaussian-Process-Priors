% collect results from both offline and online ep inference and average
% them over different data samples

clear
close all

root_path = '../../';
Results_offline = load([root_path, 'results/offline/synthetic/results_combined.mat']);
Results_online = load([root_path, 'results/online/synthetic/results_combined.mat']);

%%
undersampling_level_range = 10 : 5 : 55;
experiment_id_range = 1 : 10;

% nmse average
offline_nmse_average = zeros(1, numel(undersampling_level_range));
for i = 1 : numel(undersampling_level_range)
    undersampling_level = undersampling_level_range(i);
    offline_nmse_average(i) = nanmean(Results_offline.final_nmse{undersampling_level});
end

online_nmse_average = zeros(1, numel(undersampling_level_range));
for i = 1 : numel(undersampling_level_range)
    undersampling_level = undersampling_level_range(i);
    online_nmse_average(i) = nanmean(Results_online.final_nmse{undersampling_level});
end

% f-measure average
offline_f_measure_average = zeros(1, numel(undersampling_level_range));
for i = 1 : numel(undersampling_level_range)
    undersampling_level = undersampling_level_range(i);
    offline_f_measure_average(i) = ...
        nanmean(Results_offline.f_measure_beta_based{undersampling_level}...
        (~isnan(Results_offline.final_nmse{undersampling_level})));
end

online_f_measure_average = zeros(1, numel(undersampling_level_range));
for i = 1 : numel(undersampling_level_range)
    undersampling_level = undersampling_level_range(i);
    online_f_measure_average(i) = ...
        nanmean(Results_online.f_measure_beta_based{undersampling_level}...
        (~isnan(Results_online.final_nmse{undersampling_level})));
end

% time average
offline_time_average = zeros(1, numel(undersampling_level_range));
for i = 1 : numel(undersampling_level_range)
    undersampling_level = undersampling_level_range(i);
    offline_time_average(i) = ...
        nanmean(Results_offline.final_time{undersampling_level}...
        (~isnan(Results_offline.final_nmse{undersampling_level})));
end

online_time_average = zeros(1, numel(undersampling_level_range));
for i = 1 : numel(undersampling_level_range)
    undersampling_level = undersampling_level_range(i);
    online_time_average(i) = ...
        nanmean(Results_online.final_time{undersampling_level}...
        (~isnan(Results_online.final_nmse{undersampling_level})));
end


%%
save('synthetic_results_all.mat')
