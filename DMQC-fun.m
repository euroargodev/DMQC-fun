% DMQC-fun (Scripts to prepare and update float data for Argo DMQC.)
% 
% Here is an overview. These files are concentrated around the
% pre-calibration stage. There are more (and more updated) documentation
% in the scripts themeselves.  Note that some of them will need some
% general functions from other toolboxes, in particular EVENMAT which
% can be downloaded or cloned from
% https://github.com/evenrev1/evenmat. In order to run from your working
% directory and not from this repository, you should add the path to
% this repository's directory.
%
% init_dmqc	Name files to load/analyse, and set all paths for the 
%		reference and float data. Download of float data files from
%		the coriolis server is done here.
% load referencedata	A script to ingest, quick-check and update
%		the list of reference data for Argo DMQC (wmo_boxes.mat).
% load_floats	Build the mat-file for the float, while checking for 
%		pressure and density, existing QC-flags, etc. 
% plot_profiles	Used by load_floats to make plots if any instabilities or
%		non-monotonic increasing pressure is found.
% plot_diagnostics_ow	This is copied from the OWC toolbox and
%		modidified here, because we felt it needed some
%		improvements. It will be used instead of the
%		original/downloaded one if you run ow_calibration
%		from this local directory or make sure the path to
%		our version is first in Matlab's path list. 
%
% The following files are files in the OWC toolbox. Editing them by
% setting parameters in them are part of the recipe for calibration. In
% order not to lose your settings you can copy them to your working
% directory (or move and make symbolic links back to their original
% places). This makes updating the OWC-toolbox safer.
%
% ow_config.txt	Set the OWC paths and mapping parameters. Remember to
%		delete map- and calibration files for the float in
%		question. In my opinion the mapping parameters should
%		be separated from the path config in this file and
%		put in separate config files for each project, e.g.,
%		placed in the project folders of the float_mapped
%		folder (will fix/suggest this later).
% set_calseries	Set the calibration parameters inside this
%		script. Remember to delete calibration files
%		for the float in question. 
%
% Run the calibration by running this script from your working
% directory (but there is no need to copy it there):
%
% ow_calibration Run the (mapping) and calibration.  


% From the OWC instructions:
%
% The subdirectory /data is organised as follows (in OWC):
% /data/float_source/
% contains .mat files with the original data from the floats.
% /data/float_mapped/
% contains .mat files with the objective estimates at the float profile locations and observed θ levels.
% /data/float_calib/
% contains .mat files with the calibration constants and error estimates.
% /data/float_plots/
% contains diagnostic plots from the calibration.
% /data/constants/
% contains coastdat.mat, wmo_boxes.mat, and TypicalProfileAroundSAF.mat.
% /data/climatology/historical_ctd, /historical_bot, /argo_profiles
% are where you put your reference data. Reference data are saved as .mat files in 10×10 degree WMO boxes. Please refer to README_prepare_refdbase_ow.pdf for their data format.
 

