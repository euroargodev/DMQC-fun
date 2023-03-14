% RUN_OW_CALIBRATION is a shell to run the OW_CALIBRATION script from OWC.
% by J. Even Ø. Nilsen, Ingrid Angel, Birgit Klein, and Kjell Arne Mork.
% DMQC-fun v0.9.3, jan.even.oeie.nilsen@hi.no.

% This script will allow you to have separate ow_config.txt and
% set_calseries.m files for each float that you are analysing, and still
% run all your floats in one go. Part of the goal is also to be able to
% have automated reports for each float, incl. the mapping and
% calibration parameters, etc.
%
% This shell is necessary, because we want to keep the full list of
% floats in float_names and float_dirs, at the same time as we will run
% OW_CALIBRATION from separate diorectories for each float.
%
% DO NOT DELETE CALSERIES FILE MANUALY!  For PSAL greylisted floats,
% PREPARE_FLOATS depend on finding the calseries in set_calseries.m in
% order to keep the bad profiles in the data to be used by
% OW_CALIBRATION. 
%
% You set which floats to operate on etc. in INIT_DMQC.
%
% No editing in this file!
%
% ------------------------------------------------------------------------------


close all
init_dmqc;	% Run first to get the list and paths etc..

% Follow selection of reference data set in INIT_DMQC:
tardir=tardir(itardir); 
tartyp=tartyp(itardir);

if direction~='A'
  warning(['Direction is set to ''',direction,''' in INIT_DMQC! ',...
	   'RUN_OW_CALIBRATION is only applicable for ascending profiles. ',...
	   'Setting direction to ''A''.']); 
  direction='A';
end 

if length(tartyp)<2 % Select only one source of reference data
  switch tartyp{1}
   case {'ctd'}
    copyfile([owc_data_dir,filesep,'constants',filesep,'ctd',filesep,'wmo_boxes.mat'],[owc_data_dir,filesep,'constants',filesep,'wmo_boxes.mat'])
   case {'argo'}
    copyfile([owc_data_dir,filesep,'constants',filesep,'ctd',filesep,'wmo_boxes.mat'],[owc_data_dir,filesep,'constants',filesep,'wmo_boxes.mat'])
   otherwise
    error('Combinations other than all sources for reference data, or bottle data, not yet coded for!');
  end
end


for I=1:length(float_names)
  run([my_working_dir,'init_dmqc']);			% Run again to refresh the list.
  float_names=float_names(I);				% Use only one at a time
  float_dirs=float_dirs(I);				% Use only one at a time
  cd([my_working_dir,'DMQC',filesep,float_names{1}]);	% Go to that float's working dir
  cal_file=char(strcat(owc_data_dir,'float_calib/',float_dirs,'calseries_',float_names,'.mat'));
  delete(cal_file);					% Remove the calseries file so changes will be put to effect
  close all
  ow_calibration;					% Run on that float only (saves outcome in calib_dir)
  snippet(load_configuration('ow_config.txt'),'ow_config');		% write configuration to tex-file for report
  snippet(load(cal_file),'cal_par','',',','calseries');			% write all calibration parameters to tex-file for report
  load(cal_file); snippet(listnumstr(calseries,[],',','; '),'calseries'); % Write just the calseries to tex-file for report
end

% Make sure the pointers for the full reference data set is put back in place:
copyfile([owc_data_dir,filesep,'constants',filesep,'all',filesep,'wmo_boxes.mat'],[owc_data_dir,filesep,'constants',filesep,'wmo_boxes.mat'])
