% add EEGLAB to path
% load results from ep

time_moment = 200;

max_beta = restored_beta(:, time_moment);
active_coord_ind = find(abs(max_beta) > 300);

dipoles = ceil(active_coord_ind / 3);

sources_struct = struct;
sources_struct.posxyz = zeros(1, 3);
sources_struct.momxyz = zeros(1, 3);
sources_struct.rv = 0.1;

sources = repmat(sources_struct, 1, numel(dipoles));

for i = 1 : numel(dipoles)
    sources(i).posxyz = grid(dipoles(i), :);
    sources(i).momxyz = restored_beta((1:3) + 3 * (dipoles(i) - 1), time_moment)';
end

vol_path = '../EEGLAB/plugins/dipfit2.3/standard_BEM/standard_vol.mat';
mri_path = '../EEGLAB/plugins/dipfit2.3/standard_BEM/standard_mri.mat';

dipplot(sources, 'meshdata', vol_path, 'mri', mri_path, 'normlen', 'on', 'coordformat', 'MNI', 'projlines', 'on');
