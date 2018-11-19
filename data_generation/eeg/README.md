EEG data is generated with the help of EEGLAB and requires several steps:

1. In folder source_localisation_export run compute_leadfield_for_grid.m, which should create lead_field_matrix_and_grid.mat in that folder. This mat file will contain the lead field matrix (X) and dipole localisation grid.
2. Add folder source_localisation_export to path
3. In Matlab navigate into EEGLAB. In Maltab command line run
   ```matlab
   eeglab
   ```
4. Download ftp://sccn.ucsd.edu/pub/eeglab_dipole.set, put it into sample_data folder.
    Load dataset in EEGLAB (File -> load existing dataset).
5. In EEGLAB: Tools -> Locate dipoles using DIPFIT -> Head model -> Bounadary element model (No Co-reg)
6. In EEGLAB: Edit -> Channel locations -> Look up locs -> Use MNI
7. In Matlab navigate to EEGLAB/plugins/Fieldtrip-lite161012/private
8. In MatLab command line run
   ```matlab 
   dipfitdefs
   xgrid = eval(xgridstr);
   ygrid = eval(ygridstr);
   zgrid = eval(zgridstr);
   select = [ 1:size(EEG.icawinv,2) ];
   [ lf, grid ] = compute_leadfild_for_grid( EEG, select, xgrid, ygrid, zgrid );
   ```
9. In Matlab navigate to source_localisation_export folder and run save_event_\*.m scripts

We have 272 points(dipoles) in grid. 272 points with 3 coordinates, which are specified in grid. 
In beta:
1st dipole occupies 1, 2, 3 coordinates;
2nd dipole occupies 4, 5, 6 coordinates;
3rd dipole occupies 7, 8, 9 coordinates and so on

In total size of beta for one time moment is 816 (272 * 3)
   
