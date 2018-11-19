% collects results for all undersampling and initialisations for the
% synthetic data for online ep inference

clear
root_path = '../../';
src_path = [root_path, 'src/offline/'];
addpath(genpath(src_path)); % for nmse and f-measure functions

path=[root_path, 'results/online/synthetic/'];
% collect data
for experiment_number=1:10 
    en_str = num2str(experiment_number);
    disp(['collecting data for experiment ', en_str, '/10']);
    
    for undersampling_level=10:5:55
        ul_str = num2str(undersampling_level);
        
        SRes = load([path ,en_str,'/',ul_str,...
            '/results.mat']);
        restored_beta{experiment_number, undersampling_level} = SRes.results.restored_beta;
        time{experiment_number, undersampling_level} = SRes.results.time;
        
        clear SRes;
    end
end

% nmse and f_measure/undersampling
omega_beta_based = {};
for experiment_number = 1:10
    for undersampling_level=10:5:55
        SData = load([root_path, 'data/synthetic/', num2str(experiment_number), '/data_',num2str(undersampling_level),'_size.mat']);
        omega_beta_based{experiment_number, undersampling_level} = ...
            threshold_synthetic(restored_beta{experiment_number, undersampling_level});
        f_measure_beta_based{undersampling_level}(experiment_number) = ...
            F_measure_omega_spike(omega_beta_based{experiment_number, undersampling_level}, SData.omega);
        final_nmse{undersampling_level}(experiment_number) = ...
            NMSE(restored_beta{experiment_number, undersampling_level}, SData.beta);
    end
end

% time/undersampling
for experiment_number=1:10
    for undersampling_level=10:5:55
        final_time{undersampling_level}(experiment_number) = time{experiment_number, undersampling_level};
    end
end

save([path, 'results_combined.mat'])
