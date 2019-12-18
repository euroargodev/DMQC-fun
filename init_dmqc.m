float_names={'6902550','5904988'}; % Necessary for ow_calibration.m. Do not change name of variable!
float_names={'6902550'};
%float_names={'5904988'};
%float_names={'3902103'};
% The list of NorArgo deployment 2019:
%float_names={'6903548' '6903549' '6903550' '6903551' '6903552' '6903553' '6903554' '6903555' '6903556' '6903557' '6903558' '6903559' '6903560' '6903561' '6903562' '6903563' '6903564'};
float_names={'6903548' '6903549' '6903550' '6903551' '6903553' '6903554' '6903555' '6903557' '6903558' '6903559' '6903560' '6903561' '6903562' '6903563' '6903564'};
%float_names={'6903550'};
float_names={'6903557'}; % Central LBE with inversion.
%float_names={'6902550'};

% Root folders:
download_dir=strcat('~/Downloads/ARGO/',float_names,filesep);
rootdirin=strcat(download_dir,'Rfiles/');
rootdirout=strcat(download_dir,'Dfiles/');
% Matlabow folders:
source_dir='~/Arkiv/data/matlab_owc/float_source/';
float_dirs=repmat({'project_NorARGO/'},1,length(float_names)); % Necessary for ow_calibration.m. Do not change name of variable!
% File names:
infiles=strcat(download_dir,filesep,float_names,{'_prof.nc'});
outfiles=strcat(source_dir,float_dirs,float_names,'.mat');

MAP_P_EXCLUDE=500; % Same as in ow_config.txt but for use in our load_floats.m.


if logical(0) % ----- SWITCH ON IN ORDER TO DOWNLOAD NEW FLOAT DATA --------------------------
  for i=1:length(download_dir)
    mkdir(download_dir{i}); 
    system(['lftp -e ''lcd ',download_dir{i},' ; get ',float_names{i},'_prof.nc ; bye'' ftp.ifremer.fr/ifremer/argo/dac/coriolis/',float_names{i},'/']);
  end
end % ----------------------------------------------------------------------------------------
