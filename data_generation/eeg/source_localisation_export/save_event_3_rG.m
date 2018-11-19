% load data in eeglab

load('lead_field_matrix_and_grid.mat');

start_of_event_in_time = 0; % read event description
start_of_event_time_index = find(EEG.times == start_of_event_in_time);

epoch = 2; % read event description

y = EEG.data(:, start_of_event_time_index:end, epoch);

% navigate to saved lead field matrix

mkdir('../../../data/eeg/')
save('../../../data/eeg/data_3_event.mat', 'y', 'X', 'grid');