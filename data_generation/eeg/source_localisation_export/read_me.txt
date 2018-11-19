1. In folder source_localisation_export run compute_leadfield_for_grid.m, which should create lead_field_matrix_and_grid.mat in that folder
2. Add folder source_localisation_export to path
3. Navigate into EEGLAB. In Maltab command line run
   >> eeglab
4. Download ftp://sccn.ucsd.edu/pub/eeglab_dipole.set, put it into sample_data folder.
    Load dataset in eeglab (File -> load existing dataset).
5. Tools -> Locate dipoles using DIPFIT -> Head model -> Bounadary element model (No Co-reg)
6. Edit -> Channel locations -> Look up locs -> Use MNI
7. Navigate to EEGLAB/plugins/Fieldtrip-lite161012/private
8. In MatLab command line run
   >> dipfitdefs
   >> xgrid = eval(xgridstr);
   >> ygrid = eval(ygridstr);
   >> zgrid = eval(zgridstr);
   >> select = [ 1:size(EEG.icawinv,2) ];
   >> [ lf, grid ] = compute_leadfild_for_grid( EEG, select, xgrid, ygrid, zgrid );

9. Now you can navigate to source_localisation_export folder and run save_event_*.m scripts

We have 272 points(dipoles) in grid. 272 points with 3 coordinates, which are specified in grid. 
In beta:
1st dipole occupies 1, 2, 3 coordinates;
2nd dipole occupies 4, 5, 6 coordinates;
3rd dipole occupies 7, 8, 9 coordinates and so on

In total size of beta for one time moment is 816 (272 * 3)
   