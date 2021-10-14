% Use this shell to run the OW_CALIBRATION script.

% This script will allow you to have separate ow_config.txt and
% set_calseries.m files for each float that you are analysing, and still
% run all your floats in one go. Part of the goal is also to be able to
% have automated reports for each float, incl. the mapping and
% calibration parameters, etc.
%
% This shell is necessary, because we want to keep the full list of
% floats in float_names and float_dirs, at the same time as we will run
% OW_CALIBRATION from separate diorectories for each float.

close all
init_dmqc;	% Run first to get the list and paths etc..

for I=1:length(float_names)
  init_dmqc;						% Run again to refresh the list.
  float_names=float_names(I);				% Use only one at a time
  float_dirs=float_dirs(I);				% Use only one at a time
  cd([my_working_dir,'DMQC',filesep,float_names{1}]);	% Go to that float's working dir
  cal_file=char(strcat(owc_data_dir,'float_calib/',float_dirs,'calseries_',float_names,'.mat'));
  delete(cal_file);					% Remove the calseries file so changes will be put to effect
  close all
  ow_calibration;					% Run on that float only
  snippet(load_configuration('ow_config.txt'),'ow_config');
  snippet(load(cal_file),'cal_par','','','calseries');
  load(cal_file); [~,~,ans]=groups(calseries); snippet(ans,'calseries');
end
