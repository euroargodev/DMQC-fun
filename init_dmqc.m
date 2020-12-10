
% THE LIST OF FLOATS:
% The variable float_names is necessary for ow_calibration.m. Do not
% change name of this variable!
float_names={'6903548' '6903549' '6903550' '6903551' '6903553' '6903554' '6903555' '6903557' '6903558' '6903559' '6903560' '6903561' '6903562' '6903563' '6903564'};
corr       ={ false     false     false     true      false     false     false     false     false     false     false     false     true      false     false   };
% USE THE FOLLOWING WHILE DOING DMQC WITH OWC: 
%float_names={'6903550'}; % The one you are working on while DMQC-ing.
%float_names={'6903549','6903554','6903561','6903548'}; % The 'difficult' ones.

% DOWNLOAD directories:
download_dir=strcat('~/Downloads/ARGO/',float_names,filesep);
DIR_FTP=download_dir; % For dm_floats
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
rootdirin=strcat(download_dir,'profiles/');
rootdirout=strcat(download_dir,'Dfiles/');
% Your working directory for ARGO calibration:
my_working_dir='~/arbeid/obs/NorARGO/';

MAP_P_EXCLUDE=400; % Keep as same as in ow_config.txt but for use in our LOAD_FLOATS.


