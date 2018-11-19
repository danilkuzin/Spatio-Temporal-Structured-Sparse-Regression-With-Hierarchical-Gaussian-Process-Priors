% generates plots for quantitative comparison of results for the EEG data
% and restored vs observed EEG plots

clear
close all

root_path = '../../';
src_path = [root_path, 'src/offline/'];
addpath(genpath(src_path)); % for nmse and f-measure functions

Results_offline = load([root_path, 'results/offline/eeg/results.mat']);
Results_online = load([root_path, 'results/online/eeg/results.mat']);

Data = load([root_path, 'data/eeg/data_4_event.mat']);

%% 

% nmse over time
T = size(Data.y, 2);

nmse_over_time_offline = zeros(1, T);
for t = 1 : T
    nmse_over_time_offline(t) = norm(Data.y(:, t) - Data.X * Results_offline.results.restored_beta(:, t)) / norm(Data.y(:, t));
end

nmse_over_time_online = zeros(1, T);
for t = 1 : T
    nmse_over_time_online(t) = norm(Data.y(:, t) - Data.X * Results_online.results.restored_beta(:, t)) / norm(Data.y(:, t));
end

% plot nmse over time (around response time)
start_response_time = 163;
end_response_time = 178;

figure
plot(start_response_time:end_response_time, nmse_over_time_offline(start_response_time:end_response_time), '-', 'LineWidth', 2);
hold on
plot(start_response_time:end_response_time, nmse_over_time_online(start_response_time:end_response_time), '--', 'LineWidth', 2);
legend('offline', 'online')
saveas(gcf, 'nmse_eeg.eps', 'epsc')


%% plot EEG and EEG restored
% observed EEG
figure
imagesc(Data.y, [-0.4, 0.4])
colormap gray
colorbar
saveas(gcf, 'measured_eeg.eps', 'epsc')

figure
imagesc(Data.X * Results_offline.results.restored_beta, [-0.4, 0.4])
colormap gray
colorbar
saveas(gcf, 'restored_eeg.eps', 'epsc')


