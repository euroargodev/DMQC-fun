% The variable float_names is necessary for ow_calibration.m. Do not
% change name of this variable!
% THE LIST OF FLOATS:
float_names={'6903548' '6903549' '6903550' '6903551' '6903553' '6903554' '6903555' '6903557' '6903558' '6903559' '6903560' '6903561' '6903562' '6903563' '6903564'};
corr       ={ false     false     false     true      false     false     false     false     false     false     false     false     true      false     false   };
% USE THE FOLLOWING WHILE DOING DMQC WITH OWC: 
%float_names={'6903561'}; % The one you are working on while DMQC-ing.
%float_names={'6903549','6903554','6903561','6903548'}; % The 'difficult' ones.

% DOWNLOAD directories:
download_dir=strcat('~/Downloads/ARGO/',float_names,filesep);
download_ref_data_dir='~/Downloads/DMQC';
% MATLABOW directories:
owc_data_dir='~/Arkiv/data/matlab_owc/';
source_dir=[owc_data_dir,'float_source/'];
float_dirs=repmat({'project_NorARGO/'},1,length(float_names)); % Necessary for ow_calibration.m. Do not change name of variable!
% File names for LOAD_PROFILES:
infiles=strcat(download_dir,filesep,float_names,{'_prof.nc'}); 
outfiles=strcat(source_dir,float_dirs,float_names,'.mat'); % For load_profiles
% Directories for WRITE_DFILES:
float_dir=[owc_data_dir,'float_source',filesep,'project_NorARGO/'];
mapped_dir=[owc_data_dir,'float_mapped',filesep,'project_NorARGO/'];
calib_dir=[owc_data_dir,'float_calib',filesep,'project_NorARGO/'];
rootdirin=strcat(download_dir,'Rfiles/');
rootdirout=strcat(download_dir,'Dfiles/');

MAP_P_EXCLUDE=400; % Keep as same as in ow_config.txt but for use in our load_floats.m.

if logical(0) % ----- SWITCH ON IN ORDER TO DOWNLOAD NEW FLOAT FILES AND R-FILES ------
  for I=4%1:length(download_dir)
    mkdir(download_dir{I}); 
    system(['lftp -e ''lcd ',download_dir{I},' ; get ',float_names{I},'_prof.nc ; bye'' ftp.ifremer.fr/ifremer/argo/dac/coriolis/',float_names{I},'/']);
    mkdir(rootdirin{I}); 
    system(['lftp -e ''lcd ',rootdirin{I},' ; mget R',float_names{I},'_*.nc ; bye'' ftp.ifremer.fr/ifremer/argo/dac/coriolis/',float_names{I},'/profiles/']);
  end
end % --------------------------------------------------------------------------------
% R-files must be downloaded at the same time, or else there might be
% mismatch between number of cal parameters and R-files.
