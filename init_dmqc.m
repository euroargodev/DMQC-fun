% Variable float_names is necessary for ow_calibration.m. Do not change name of variable!
%float_names={'6902550','5904988'}; 
%float_names={'6902550'};
%float_names={'5904988'};
%float_names={'3902103'};
% The list of NorArgo deployment 2019:
%float_names={'6903548' '6903549' '6903550' '6903551' '6903552' '6903553' '6903554' '6903555' '6903556' '6903557' '6903558' '6903559' '6903560' '6903561' '6903562' '6903563' '6903564'};
%float_names={'6903550'};
%float_names={'6903557'}; % Central LBE with inversion.
%float_names={'6902550'};
% NorARGO
float_names={'6903548' '6903549' '6903550' '6903551' '6903553' '6903554' '6903555' '6903557' '6903558' '6903559' '6903560' '6903561' '6903562' '6903563' '6903564'};
%float_names={'6903560'};
%float_names={'6903548'}; % The one you are working on 

% Root folders:
download_dir=strcat('~/Downloads/ARGO/',float_names,filesep);
% Matlabow folders:
owc_data_dir='~/Arkiv/data/matlab_owc/';
source_dir=[owc_data_dir,'float_source/'];
float_dirs=repmat({'project_NorARGO/'},1,length(float_names)); % Necessary for ow_calibration.m. Do not change name of variable!
% Matlabow folders for WRITE_DFILES:
float_dir=[owc_data_dir,'float_source',filesep,'project_NorARGO/'];
mapped_dir=[owc_data_dir,'float_mapped',filesep,'project_NorARGO/'];
calib_dir=[owc_data_dir,'float_calib',filesep,'project_NorARGO/'];
rootdirin=strcat(download_dir,filesep,'Rfiles/');
rootdirout=strcat(download_dir,filesep,'Dfiles/');
% File names for LOAD_PROFILES:
infiles=strcat(download_dir,filesep,float_names,{'_prof.nc'}); 
outfiles=strcat(source_dir,float_dirs,float_names,'.mat'); % For load_profiles

MAP_P_EXCLUDE=400; % Same as in ow_config.txt but for use in our load_floats.m.


if logical(0) % ----- SWITCH ON IN ORDER TO DOWNLOAD NEW FLOAT DATA --------------------------
  for i=1:length(download_dir)
    mkdir(download_dir{i}); 
    system(['lftp -e ''lcd ',download_dir{i},' ; get ',float_names{i},'_prof.nc ; bye'' ftp.ifremer.fr/ifremer/argo/dac/coriolis/',float_names{i},'/']);
  end
end % ----------------------------------------------------------------------------------------
