% RUN_OW_CALIBRATION is a shell to run the OW_CALIBRATION script from
% the OWC toolbox on several floats.
%
% DMQC-fun v0.9.
% J. Even Ø. Nilsen, Ingrid M. Angel-Benavides, Birgit Klein, Malgorzata Merchel, and Kjell Arne Mork.
% Last updated: Wed May 24 13:40:44 2023 by jan.even.oeie.nilsen@hi.no

% This script will allow you to have separate ow_config.txt and
% set_calseries.m files for each float that you are analysing, and still
% run all your floats in one go. Part of the goal is also to be able to
% have automated reports for each float, incl. the mapping and
% calibration parameters, etc.
%
% This shell is necessary, because we want to keep the full list of
% floats in float_names and float_dirs, and run OW_CALIBRATION from
% separate diorectories for each float.
%
% Delete float_mapped files manually if you want to redo the OWC
% mapping. The float_calib files are automatically replaced here.
% 
% DO NOT DELETE CALSERIES FILE MANUALY!  For PSAL greylisted floats,
% PREPARE_FLOATS depend on finding the calseries in set_calseries.m in
% order to keep the bad profiles in the data to be used by
% OW_CALIBRATION, and you can also set the clustering in
% PREPARE_FLOATS to follow the calseries.
%
% You set which floats to operate on etc. in INIT_DMQC.
%
% List of actions in this script:
% - runs OW_CALIBRATION from the OWC toolbox on each float.
% - writes snippets about OWC-configuration for the report.
% - makes the whole section about OWC for the report.
% - erases the main section from PAIR_FLOATS_WITH_REFERENCEDATA.
% - if OWC fails:
%	- only makes the map figure from OWC for the report.
%	- makes dummy float_calib files for WRITE_D.
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
  tic; 
  run([my_working_dir,'init_dmqc']);			% Run again to refresh the list.
  float_names=float_names(I);				% FOR OW_CALIBRATION, use only one at a time
  float_dirs=float_dirs(I);				% FOR OW_CALIBRATION, use only one at a time
  cal_action=cal_action(I);				% local, use only one at a time
  outfiles=outfiles(I);					% local, use only one at a time
  cd([my_working_dir,'DMQC',filesep,float_names{1}]);	% Go to that float's working dir
  calseries_file = char(strcat(calib_dir,'calseries_',float_names,'.mat'));
  cal_file       = char(strcat(calib_dir,'cal_',      float_names,'.mat'));
  if exist(calseries_file,'file'), delete(calseries_file); end	% Remove the calseries file so changes will be put to effect
  if exist(cal_file,'file'),       delete(cal_file);       end	% Remove the cal file so changes will be put to effect
  close all
  try
    ow_calibration;								% Run on that float only (saves outcome in calib_dir)
    snippet(load_configuration('ow_config.txt'),'ow_config');			% write configuration to tex-file for report
    snippet(load(calseries_file),'cal_par','',',','calseries');			% write all calibration parameters to tex-file for report
    load(calseries_file); snippet(listnumstr(calseries,[],',','; '),'calseries'); % Write just the calseries to tex-file for report
    ct_run_ow_calibration=toc/60
    
    % Make section for report:
    fid=fopen('owc.tex','w');
    fprintf(fid,'%s\n','\subsubsection{Salinity calibration with OWC}');
    fprintf(fid,'%s\n','Figure~\ref{trajectory} shows the positions of reference data actually used in mapping.');
    fprintf(fid,'%s\n','For floats with previously assigned uncorrectable profiles, these profiles are still included in the current OWC analysis');
    fprintf(fid,'%s\n','in order to avoid gaps in the presentation of the mapping and calibration (esp.\ Figure~\ref{SalWithErrors}). ');
    fprintf(fid,'%s\n','');
    fprintf(fid,'%s\n','The following are the mapping configuration parameters set in \texttt{ow\_config.txt} file of the OWC toolbox with the parameters set for the final correction:  ');
    fprintf(fid,'%s\n','\verbatiminput{ow_config}');
    fprintf(fid,'%s\n','');
    fprintf(fid,'%s\n','The calseries parameters are set in \texttt{set\_calseries.m} file as follows: ');
    fprintf(fid,'%s\n','\verbatiminput{cal_par}');
    fprintf(fid,'%s\n','');
    fprintf(fid,'%s\n','Figures \ref{trajectory} through \ref{Salinity_OWlevels} show the results of the comparison and correction of the salinity data.');
    fprintf(fid,'%s\n','Notes made about this float during the different rounds of DMQC, can be found in Appendix~A.');
    fprintf(fid,'%s\n','');
    fprintf(fid,'%s\n','\begin{figure}[ph]');
    fprintf(fid,'%s\n','  \centering    ');
    fprintf(fid,'%s\n','  \includegraphics[width=\textwidth,natwidth=450,natheight=330]{\floatplots\WMOnum_1.eps}');
    fprintf(fid,'%s\n','  \caption{Float \WMOnum. Location of the float profiles (red line with');
    fprintf(fid,'%s\n','    black numbers) and the reference data selected for mapping (blue');
    fprintf(fid,'%s\n','    dots).');
    fprintf(fid,'%s\n','    Note that missing cycles are not included into OWC, so profile numbers may');
    fprintf(fid,'%s\n','    not match cycle numbers.}');
    fprintf(fid,'%s\n','  \label{trajectory}');
    fprintf(fid,'%s\n','\end{figure}');
    fprintf(fid,'%s\n','\begin{figure}[ph]');
    fprintf(fid,'%s\n','    \centering    ');
    fprintf(fid,'%s\n','    \includegraphics[width=\textwidth,natwidth=1400,natheight=820]{\floatplots\WMOnum_3}');
    fprintf(fid,'%s\n','    \caption{Float \WMOnum. Evolution of the suggested adjustment with');
    fprintf(fid,'%s\n','      time. The top panel plots the potential conductivity multiplicative');
    fprintf(fid,'%s\n','      adjustment. The bottom panel plots the equivalent salinity additive');
    fprintf(fid,'%s\n','      adjustment. The red line denotes one-to-one profile fit that uses the');
    fprintf(fid,'%s\n','      vertically weighted mean of each profile. The red line can be used to');
    fprintf(fid,'%s\n','      check for anomalous profiles relative to the optimal fit.');
    fprintf(fid,'%s\n','      Note that missing cycles are not included into OWC, so profile numbers do');
    fprintf(fid,'%s\n','      not necessarily match cycle numbers.} ');
    fprintf(fid,'%s\n','    \label{SalWithErrors}');
    fprintf(fid,'%s\n','\end{figure}');
    fprintf(fid,'%s\n','\begin{figure}[p]');
    fprintf(fid,'%s\n','    \centering    ');
    fprintf(fid,'%s\n','    \includegraphics[width=\textwidth,natwidth=1000,natheight=900]{\floatplots\WMOnum_2}');
    fprintf(fid,'%s\n','    \caption{Float \WMOnum. The original float salinity and the');
    fprintf(fid,'%s\n','      objectively estimated reference salinity at the 10 float theta');
    fprintf(fid,'%s\n','      levels that are used in calibration as errorbars. Lower panel is a');
    fprintf(fid,'%s\n','      zoom to the latter.');
    fprintf(fid,'%s\n','      Note that missing cycles are not included into OWC, so profile');
    fprintf(fid,'%s\n','      numbers in legend do not necessarily match cycle numbers.}');
    fprintf(fid,'%s\n','    \label{uncalibVsSalinity}');
    fprintf(fid,'%s\n','\end{figure}');
    fprintf(fid,'%s\n','\begin{figure}[p]');
    fprintf(fid,'%s\n','    \centering    ');
    fprintf(fid,'%s\n','    \includegraphics[width=\textwidth,natwidth=1000,natheight=900]{\floatplots\WMOnum_4}');
    fprintf(fid,'%s\n','    \caption{Float \WMOnum. Plots of calibrated float salinity and the');
    fprintf(fid,'%s\n','      objectively estimated reference salinity at the 10 float theta');
    fprintf(fid,'%s\n','      levels that are used in calibration as errorbars. Lower panel is a');
    fprintf(fid,'%s\n','      zoom to the latter.');
    fprintf(fid,'%s\n','      Note that missing cycles are not included into OWC, so profile');
    fprintf(fid,'%s\n','      numbers in legend do not necessarily match cycle numbers.}');
    fprintf(fid,'%s\n','    \label{CalibVsSalinity}');
    fprintf(fid,'%s\n','\end{figure}');
    fprintf(fid,'%s\n','\begin{figure}[p]');
    fprintf(fid,'%s\n','    \centering     ');
    fprintf(fid,'%s\n','    \includegraphics[width=0.8\textwidth,natwidth=500,natheight=700]{\floatplots\WMOnum_5}');
    fprintf(fid,'%s\n','    \caption{Float \WMOnum. Salinity anomaly on theta levels.');
    fprintf(fid,'%s\n','      Note that missing cycles are not included into OWC, so profile numbers do');
    fprintf(fid,'%s\n','      not necessarily match cycle numbers.}');
    fprintf(fid,'%s\n','    \label{SalAnomOnTheta}');
    fprintf(fid,'%s\n','\end{figure}');
    fprintf(fid,'%s\n','\begin{figure}[p]');
    fprintf(fid,'%s\n','    \centering    ');
    fprintf(fid,'%s\n','    \includegraphics[width=0.8\textwidth,natwidth=500,natheight=700]{\floatplots\WMOnum_7}');
    fprintf(fid,'%s\n','    \caption{Float \WMOnum.  C alibrated salinity anomaly on theta levels.');
    fprintf(fid,'%s\n','      Note that missing cycles are not included into OWC, so profile numbers do');
    fprintf(fid,'%s\n','      not necessarily match cycle numbers.}');
    fprintf(fid,'%s\n','    \label{CalibSalAnomOnTheta}');
    fprintf(fid,'%s\n','\end{figure}');
    fprintf(fid,'%s\n','\begin{figure}[p]');
    fprintf(fid,'%s\n','    \centering    ');
    fprintf(fid,'%s\n','    \includegraphics[width=\textwidth,natwidth=1350,natheight=820]{\floatplots\WMOnum_6}');
    fprintf(fid,'%s\n','    \caption{Float \WMOnum. Plots of the evolution of salinity with time');
    fprintf(fid,'%s\n','      along with selected theta levels with minimum salinity variance.');
    fprintf(fid,'%s\n','      Note that missing cycles are not included into OWC, so profile numbers do');
    fprintf(fid,'%s\n','      not necessarily match cycle numbers.} ');
    fprintf(fid,'%s\n','    \label{SalErrOnTheta}');
    fprintf(fid,'%s\n','\end{figure}');
    fprintf(fid,'%s\n','\begin{figure}[p]');
    fprintf(fid,'%s\n','  \centerline{\includegraphics[width=1.2\textwidth,natwidth=1700,natheight=1000]{\floatplots\WMOnum_8}}');
    fprintf(fid,'%s\n','  \caption{Float \WMOnum. Plots include the theta levels chosen for');
    fprintf(fid,'%s\n','    calibration: Top left: Salinity variance at theta levels. Top right:');
    fprintf(fid,'%s\n','    T/S diagram of all profiles of Argo float. Bottom left: potential');
    fprintf(fid,'%s\n','    temperature plotted against pressure. Bottom right: salinity plotted');
    fprintf(fid,'%s\n','    against pressure.} ');
    fprintf(fid,'%s\n','  \label{Salinity_OWlevels}');
    fprintf(fid,'%s\n','\end{figure}');
    fprintf(fid,'%s\n','');
    fclose(fid); 
    
    % If OWC is sucessful, the coinciding comparisons are superfluous: 
    if ~exist('coincidence.tex','file'), fid=fopen('coincidence.tex','w'); fclose(fid); end
    % However, we do not empty the coincidence_appendix.tex if it has been made.
   
  catch % if OWC is NOT sucessful: 
    disp([float_names{1},': ow_calibration failed! No section for the report (only map-figure in owc.tex).']);
    % Make file for report with the map figure only (no section title):
    fid=fopen('owc.tex','w');
    fprintf(fid,'%s\n','Figure~\ref{trajectory} shows the trajectory and where there are nearby reference data.');
    fprintf(fid,'%s\n','\begin{figure}[ph]');
    fprintf(fid,'%s\n','  \centering    ');
    fprintf(fid,'%s\n','  \includegraphics[width=0.8\textwidth,natwidth=450,natheight=330]{\floatplots\WMOnum_1.eps}');
    fprintf(fid,'%s\n','  \caption{Float \WMOnum. Location of the float profiles (red line with black numbers) and the reference data selected by the attempted OWC-mapping (blue dots).');
    fprintf(fid,'%s\n','    Note that missing cycles are not included into OWC, so profile numbers may not match cycle numbers.}');
    fprintf(fid,'%s\n','  \label{trajectory}');
    fprintf(fid,'%s\n','\end{figure}');
    fclose(fid); 
    % % Make empty file for report:
    % fid=fopen('owc.tex','w'); fclose(fid); 
    % Make passive calseries object for write_D:
    load(outfiles{1},'SAL','PTMP');
    [m,n]=size(PTMP);
    calseries=ones(1,size(PTMP,2));
    save(calseries_file,'calseries'); 
    % Make passive cal_ file objects for write_D:
    pcond_factor=1;
    pcond_factor_err=1e-6;
    cal_COND = sw_c3515*sw_cndr(SAL,PTMP,0);
    cal_COND_err = pcond_factor_err.*cal_COND;
    cal_SAL=SAL;
    cal_SAL_err = abs(cal_SAL-sw_salt((cal_COND+cal_COND_err)/sw_c3515,PTMP,0));
    sta_SAL = nan(m,n);
    sta_SAL_err = nan(m,n);
    save(cal_file,'cal_COND*','cal_SAL*','pcond_*','sta_SAL*'); 
  end

end

% Make sure the pointers for the full reference data set is put back in place:
copyfile([owc_data_dir,filesep,'constants',filesep,'all',filesep,'wmo_boxes.mat'],[owc_data_dir,filesep,'constants',filesep,'wmo_boxes.mat'])
