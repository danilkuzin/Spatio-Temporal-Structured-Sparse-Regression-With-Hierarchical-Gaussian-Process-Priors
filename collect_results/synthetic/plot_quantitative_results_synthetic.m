% generates plots for quantitative comparison of results for the synthetic 
% data

clear;
close all;

load('synthetic_results_all.mat');

% plot logarithmic scale nmse
figure
semilogy(undersampling_level_range, offline_nmse_average)
hold on
semilogy(undersampling_level_range, online_nmse_average)
legend('offline', 'online')
title('nmse logarithmic scale')
saveas(gcf, 'nmse_synthetic.eps', 'epsc')

% plot f measure
figure
plot(undersampling_level_range, offline_f_measure_average)
hold on
plot(undersampling_level_range, online_f_measure_average)
legend('offline', 'online')
title('f measure')
saveas(gcf, 'f_measure_synthetic.eps', 'epsc')

% plot time
figure
plot(undersampling_level_range, offline_time_average)
hold on
plot(undersampling_level_range, online_time_average)
legend('offline', 'online')
title('time')
saveas(gcf, 'time_synthetic.eps', 'epsc')
