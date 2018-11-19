% generates plots and video of restored dipoles for the EEG data. NOTE that
% video generation would generate a lot of popping up figures

clear
close all

%% load results
root_path = '../../';
addpath(genpath([root_path, 'data_generation/eeg/EEGLAB']))
load([root_path, 'results/offline/eeg/results.mat'])
load([root_path, 'data/eeg/data_4_event.mat'])

%% snapshots at specific time moments
selected_time_moments = [1, 170]; % no response yet, response time
for time_moment = selected_time_moments
    max_beta = results.restored_beta(:, time_moment);
    active_coord_ind = find(abs(max_beta) > 400);

    dipoles = ceil(active_coord_ind / 3);

    sources_struct = struct;
    sources_struct.posxyz = zeros(1, 3);
    sources_struct.momxyz = zeros(1, 3);
    sources_struct.rv = 0.1;

    sources = repmat(sources_struct, 1, numel(dipoles));

    for i = 1 : numel(dipoles)
        sources(i).posxyz = grid(dipoles(i), :);
        sources(i).momxyz = results.restored_beta((1:3) + 3 * (dipoles(i) - 1), time_moment)';
    end

    vol_path = [root_path, 'data_generation/eeg/EEGLAB/plugins/dipfit2.3/standard_BEM/standard_vol.mat'];
    mri_path = [root_path, 'data_generation/eeg/EEGLAB/plugins/dipfit2.3/standard_BEM/standard_mri.mat'];

    dipplot(sources, 'meshdata', vol_path, 'mri', mri_path, 'normlen', 'on', 'coordformat', 'MNI', 'projlines', 'on', 'view', [1 -1 1]);

    drawnow;

    F = getframe(gca);

    dipole_image = F.cdata(:, 163:end, :);

    figure;
    imshow(dipole_image);

    saveas(gcf, ['dipole_time_moment_', num2str(time_moment), '.eps'], 'epsc');
end


%% video dipole and eeg input

videoWriter = VideoWriter('eeg_source_results_event_4_and_eeg_input.mp4', 'MPEG-4');
videoWriter.FrameRate = 30;
open(videoWriter);

[K, T] = size(y);

for time_moment = 1 : T
    max_beta = results.restored_beta(:, time_moment);
    active_coord_ind = find(abs(max_beta) > 400);

    dipoles = ceil(active_coord_ind / 3);

    sources_struct = struct;
    sources_struct.posxyz = zeros(1, 3);
    sources_struct.momxyz = zeros(1, 3);
    sources_struct.rv = 0.1;

    sources = repmat(sources_struct, 1, numel(dipoles));

    for i = 1 : numel(dipoles)
        sources(i).posxyz = grid(dipoles(i), :);
        sources(i).momxyz = results.restored_beta((1:3) + 3 * (dipoles(i) - 1), time_moment)';
    end

    vol_path = [root_path, 'data_generation/eeg/EEGLAB/plugins/dipfit2.3/standard_BEM/standard_vol.mat'];
    mri_path = [root_path, 'data_generation/eeg/EEGLAB/plugins/dipfit2.3/standard_BEM/standard_mri.mat'];

    dipplot(sources, 'meshdata', vol_path, 'mri', mri_path, 'normlen', 'on', 'coordformat', 'MNI', 'projlines', 'on', 'view', [1 -1 1]);

    
    drawnow;
    
    F = getframe(gca);
    
    dipole_image = F.cdata(:, 163:end, :);
    
    FigHandle = figure('Position', [100, 100, 900, 1000]);
    subplot('position', [0 0.4 1 0.6])
    imshow(dipole_image);
    title('Dipoles');
    subplot('position', [0.05 0.05 0.9 0.35])
    ax = gca;
    ax.FontSize = 8;
    ylim([0 14]);
    xlim([0 T]);
    set(ax,'ytick',[])
    xlabel('Time');
    hold on
    j = K : -1 : 1;
    for i = K : -1 : 1
        plot(y(i, 1:time_moment) + j(i) * 0.2, 'b');
    end
    
    
    f_out = getframe(FigHandle);
    writeVideo(videoWriter, f_out);
    
    close all;
    
end


close(videoWriter);





