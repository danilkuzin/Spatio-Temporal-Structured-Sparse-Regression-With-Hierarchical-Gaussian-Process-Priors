# Spatio-Temporal-Structured-Sparse-Regression-With-Hierarchical-Gaussian-Process-Priors

The code to reproduce results from the paper ["Spatio-Temporal Structured Sparse Regression With Hierarchical Gaussian Process Priors"](https://arxiv.org/abs/1807.05561) by Danil Kuzin, Olga Isupova, and Lyudmila Mihaylova in IEEE Transactions on Signal Processing (2018) (offline and online EP inference for the proposed model for the synthetic and EEG data).

The main algorithms for offline and online EP inference are implemented in src folder.

## Data generation
Folder data_generation contains scripts to generate the synthetic and EEG data used in the paper to a newly generated data folder. 

* The synthetic data generates directly from data_generation/synthetic/generate_data.m

* The EEG data requires EEGLAB (added as submodule to this repository) and some manual steps described in data_generation/eeg/README.md

## Experiment run
Folder experiments contains scripts to run offline and online EP inference for the proposed model from the paper on the synthetic and EEG data. They store results into a newly generated results folder.

Scripts run_offline.m and run_online.m are functions to run offline and online EP inference for the proposed model respectively on any given data. 

## Plotting
Folder collect_results contains scripts to generate various plots on results stored in the results folder. These are plots from the paper and a video of active dipoles found in the EEG data.
