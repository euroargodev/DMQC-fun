float_names={'6902550','5904988'};			% Necessary for ow_calibration.m. Do not change!
float_names={'6902550'};
%float_names={'5904988'};
%float_names={'3902103'};

% The list of NorArgo deployment 2019:
float_names={'6903548' '6903549' '6903550' '6903551' '6903552' '6903553' '6903554' '6903555' '6903556' '6903557' '6903558' '6903559' '6903560' '6903561' '6903562' '6903563' '6903564'};
float_names={'6903548' '6903549' '6903550' '6903551' '6903553' '6903554' '6903555' '6903557' '6903558' '6903559' '6903560' '6903561' '6903562' '6903563' '6903564'};
%float_names={'6903550'};
%float_names={'6903557'}; % Central LBE with inversion.
%float_names={'6902550'};

% Root folders:
download_dir=strcat('~/Downloads/ARGO/',float_names,filesep);
rootdirin=strcat(download_dir,'Rfiles/');
rootdirout=strcat(download_dir,'Dfiles/');
% Matlabow folders:
source_dir='~/Arkiv/data/matlabow/float_source/';
float_dirs=repmat({'project_NorARGO/'},1,length(float_names)); % Necessary for ow_calibration.m. Do not change!
% Files:
infiles=strcat(download_dir,filesep,float_names,{'_prof.nc'});
outfiles=strcat(source_dir,float_dirs,float_names,'.mat');

MAP_P_EXCLUDE=200;


if logical(0) % ----- SWITCH ON IN ORDER TO DOWNLOAD NEW FLOAT DATA --------------------------
  for i=1:length(download_dir)
    mkdir(download_dir{i}); 
    system(['lftp -e ''lcd ',download_dir{i},' ; get ',float_names{i},'_prof.nc ; bye'' ftp.ifremer.fr/ifremer/argo/dac/coriolis/',float_names{i},'/']);
  end
end % ----------------------------------------------------------------------------------------




% Here is an overview of the files Even has made for a start. These
% files are concentrated around the pre-calibration stage. There are
% more documentation in the scripts themeselves.
%
% Note that some of them will need small functions from other toolboxes,
% in particular EVENMAT which can be downloaded at
% https://github.com/evenrev1/evenmat .
%
% init_dmqc	Name files to load/analyse, and set all paths for the 
%		DMQC. Download NetCDF file from coriolis server is
%		done here
% load referencedata	A script to ingest, quick-check and update
%		list of reference data for Argo DMQC.
% load_floats	Build the mat-file for the float, while checking for 
%		pressure and density, and existing QC-flags. 
% plot_profiles	Used by load_floats to make plots if any instabilities or
%		non-monotonic increasing pressure found.
% plot_diagnostics_ow	This is copied from the OWC toolbox and
%		modidified here, because we felt it needed some
%		improvements. It will be used instead of the
%		original/downloaded one if you run ow_calibration
%		from this local directory or make sure the path to
%		our version is first in Matlab's path list. 
%
% The following files are files in the OWC toolbox. I just mention these
% for completeness, as editing them by setting parameters in them are
% part of the recipe for calibration:
%
% ow_config.txt	Set the OWC paths and mapping parameters. Remember to
%		delete map- and calibration files for the float in
%		question. 
% set_calseries	Set the calibration parameters inside this
%		script. Remember to delete map- and calibration files
%		for the float in question. 
% ow_calibration Run the (mapping) and calibration.  

