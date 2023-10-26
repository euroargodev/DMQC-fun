% PREPARE_FLOATS loads float nc-files from your download folder, checks
% and prepares float data for ow_calibration, as well as adding flags
% and other DMQC parameters.
%
% DMQC-fun v0.9.
% J. Even Ø. Nilsen, Ingrid M. Angel-Benavides, Birgit Klein, Malgorzata Merchel, and Kjell Arne Mork.
% Last updated: Thu Oct 26 15:14:12 2023 by jan.even.oeie.nilsen@hi.no
%
% You will likely run this at least twice when new profiles come in.
% You set which floats to operate on etc. in INIT_DMQC. You also select
% direction in INIT_DMQC, in order to inspect both ascending and
% descending profiles.
%
% This script basicly follows the Argo QC manual, and does the following:
% - Ingests individual R-files
% - Loads QC flags from previous DMQC sessions (see below about flagging)
% - Fills scientific calib information
% - Delayed-mode procedures for coordinates and pressure
% - Pressure adjustments for APEX floats
% - Delayed-mode procedures for temperature and salinity:
%	Automated tests prior to visual control
%	Facilitates visual control of all profiles and flags
% - Makes several graphics necessary for DMQC
% - Applies flags to the data (in matlab)
% - Corrects deep argo/arvor floats pressure dependent conductivity bias
% - Saves results for both OWC and WRITE_D
% - Creates several parts of the report:
%	snippets of text, e.g., for the technical table
%	altimetry section
%	reference data appendix
%	visual test example figures
%	Overall comparison with the reference data figures
%	Sections about the APEX and deep corrections
%	The flagging summary table
%
% No editing in this file!
%
% ----- On flagging throughout PREPARE_FLOATS and WRITE_D: -------------
%
% The flag matrices at work here are: 
%
% *qco = old flags collected from the R-files, and for previously
%        DMQC'ed profiles filled with flags from previously read
%        R-files, thus *qco is always RTQC flags. 
%
% *qcn = new flags collected from previously saved DMQC or just
%        initiated as blank for new profiles. DMQC flags of new (and
%        possibly old) profiles are added to this if changed in a
%        session. New flags can also be changed in a later session. 
%
% *qct = temporary flags based on *qco flags, to which flags from
%        automatic and visual DMQC is added, before being assigned
%        to *qcn if they are different from flags in *qco and
%        *qcn. Flags in *qct are only used by the automatic and
%        visual controls within this script. 
% 
% Remember that you can choose to rerun DMQC on any profiles by setting
% the 'checked' parameter in INIT_DMQC back to whichever profile you
% want to start DMQC from, and stop during visual control whenever you
% like. Note that you will then start from scratch on those profiles and
% see only RTQC flags and the automatic checks. Flags from previous DMQC
% will not show nor carry through to *qcn for those particular profiles.
%
% The flags in *qco and *qcn will be kept apart at all stages in order
% to keep the statistics of RTQC and DMQC for the whole life of the
% float. They are saved in a file in the float_working_dir (a directory
% your system should keep backups of) so that one can separate these
% through several instances of DMQC. This is necessary for good
% statistics in the report, as D-files do not discern between RTQC and
% DMQC flags (Note that only DMQC flags of type '4' and '1' enter
% D-files).
% 
% However, for floats DMQC'ed before using this version of DMQC-fun the
% DMQC flags will have been written to the D-files back then, and thus
% RTQC and DMQC flags of previously inspected profiles will look the
% same. Thus, for plot marks and statistics in PREPARE_FLOATS older DMQC
% flags will be marked and counted as RTQC flags.
%
% ----- On columns, profiles, and cycles -------------------------------
% 
% Each column of the matrices used herein corresponds to a file in the
% Rfiles list, which is not (necessarily) the same as cycle numbers. But
% since the Rfiles list is ported to WRITE_D and used there, the output
% of columns from the matrices will be put in the correct D-files. But
% do take care to use the variable CYCLE_NUMBER throughout when needed
% in plots and text!
%
% ----------------------------------------------------------------------

clear all; close all
init_dmqc; % Paths and filenames.
grey=[.5 .5 .5]; % Colour setting for maps.

% Follow selection of reference data set in INIT_DMQC:
tardir=tardir(itardir); 
tartyp=tartyp(itardir);

for I=1:length(download_dir)	% Loop floats

  close all
  clear hP hPo hPn LO LA scientific_calib
  
  % ----- INITIAL INFO EXTRACTION: --------------------------------------

  % Specifics for ascending and descending profiles:
  switch direction
   case 'A',
    disp(['prepare_floats: Processing ascending profiles from float ',float_names{I},'.']);
    float_working_dir=[my_working_dir,'DMQC',filesep,float_names{I}];
    checked_name='checked';
   case 'D'
    disp(['prepare_floats: Processing DESCENDING profiles from float ',float_names{I},'.']);
    float_working_dir=[my_working_dir,'DMQC',filesep,float_names{I},filesep,'descending'];
    if ~exist(float_working_dir,'dir'), mkdir(float_working_dir); end
    outfiles{I}=[outfiles{I}(1:end-4),'D',outfiles{I}(end-3:end)];
    checked{I}=Dchecked{I};
    checked_name='Dchecked';
    % For descending profiles all snippets, tex files (fprintf) and the
    % cluster-data mat-files will go in the descending subdirectory. The
    % figures (print) and final saved mat-file with results will all be
    % marked D like the input.
  end
  
  % Go to that float's correct working dir:
  cd(float_working_dir);

  % Extract the MAP_P_EXCLUDE from the float-specific setup file:
  aa=strip(readlines([my_working_dir,'DMQC',filesep,float_names{I},filesep,'ow_config.txt']));	% Read local ow_config.txt
  find(contains(aa,'MAP_P_EXCLUDE=') & ~startsWith(aa,'%')); % Scan string 
  eval([char(aa(ans(end))),';']);			% Get the last definition of MAP_P_EXCLUDE

  % Generate snippets of text for the report (from the high level files):
  snippet(refdir);
  ncread(metafiles{I},'DATA_CENTRE'); 
  switch ans'
   case 'IF', 'Coriolis';
  end
  snippet(ans,'DAC','_deblank');
  ncread(metafiles{I},'FLOAT_SERIAL_NO'); snippet(ans','float-serial-no','_deblank');
  platyp=ncread(metafiles{I},'PLATFORM_TYPE'); platyp=snippet(platyp','platform-type','_deblank');
  ncread(metafiles{I},'TRANS_SYSTEM'); snippet(ans(:,1)','transmission-system','_deblank');
  sensor=ncread(metafiles{I},'SENSOR')';
  sensormodel=ncread(metafiles{I},'SENSOR_MODEL')'; 
  sensorsnr=ncread(metafiles{I},'SENSOR_SERIAL_NO')'; 
  contains(string(sensor),'CTD_TEMP')|contains(string(sensor),'CTD_CNDC'); ctdi=find(ans); othi=find(~ans);
  ctdmodel=snippet(sensormodel(ctdi(1),:),'ctdmodel','_deblank');
  snippet(sensorsnr(ctdi(1),:),'ctdsensorsnr','_deblank');
  snippet(sensor(othi,:),'othersensors','_deblank');
  snippet(sensormodel(othi,:),'othersensormodels','_deblank');
  snippet(sensorsnr(othi,:),'othersensorsnr','_deblank');
  ncread(metafiles{I},'LAUNCH_DATE'); snippet([ans(7:8)','/',ans(5:6)','/',ans(1:4)'],'depdate');
  ncread(metafiles{I},'LAUNCH_LATITUDE'); snippet(ans,'deplat');
  ncread(metafiles{I},'LAUNCH_LONGITUDE'); snippet(ans,'deplon');
  ncread(trajfiles{I},'REPRESENTATIVE_PARK_PRESSURE'); snippet([int2str(round(nanmedian(ans)/100)*100),' m'],'park-depth');
  %ncread(trajfiles{I},'REPRESENTATIVE_PARK_PRESSURE'); snippet([int2str(round(max(ans(end))/100)*100),' m'],'park-depth');
  %ncread(proffiles{I},'PRES'); snippet([int2str(round(max(ans(:,end))/100)*100),' m'],'profile-depth');
  ncread(proffiles{I},'PRES'); 
  profdepth=round(nanmax(ans,[],'all')/100)*100; snippet([int2str(profdepth),' m'],'profile-depth');
  ncread(trajfiles{I},'JULD_ASCENT_END'); t=ans; snippet([int2str(round(nanmedian(diff(t,1,1)))),' days'],'cycle-time');
  %ncread(proffiles{I},'JULD'); t=ans; snippet([int2str(round(nanmean(diff(t,1,1)))),' days'],'cycle-time');
  ncread(metafiles{I},'DEPLOYMENT_PLATFORM'); snippet(ans(:,1)','ship','_deblank');
  ncread(metafiles{I},'PI_NAME'); snippet(ans(:,1)','PI','deblank');
  ncread(metafiles{I},'END_MISSION_STATUS'); 
  switch ans'
   case ' ', 'Active';
   case 'T', 'No more transmission received';
   case 'R', 'Retrieved';
  end
  snippet(ans,'float-status','_deblank');
  snippet([num2str(diff(mima(t))/365,'%4.1f'),' yrs'],'age');
  inCYCLE_NUMBER=ncread(infiles{I},'CYCLE_NUMBER')'; snippet(max(inCYCLE_NUMBER),'last-cycle');
  snippet(zipnumstr(setdiff(1:max(inCYCLE_NUMBER),inCYCLE_NUMBER)),'missing-cycles');
  % Grey list:
  [g1,g2,g3,g4,g5,g6,g7]=textread([download_parent_dir,'coriolis_greylist.csv'],'%s%s%s%s%s%s%s','delimiter',',');
  ii=find(contains(g1,float_names{I}));
  glist=''; for i=1:length(ii), glist=[glist,g2{ii(i)},',',g5{ii(i)},',',g6{ii(i)},';  ']; end
  snippet(glist,'grey-list','_');
  % The set minimum depth in OWC:
  snippet(MAP_P_EXCLUDE);
  % Transfer the operator info to the report:
  snippet(dmqc_operator(1).name,'dmqc_operator_name');
  snippet(dmqc_operator(1).orcid,'dmqc_operator_orcid');
  snippet(dmqc_operator(1).institution,'dmqc_operator_institution');
  snippet(dmqc_operator(1).address,'dmqc_operator_address');
  % Make some necessary definitions for LaTeX report: 
  fid=fopen('wmonumdef.tex','w'); fprintf(fid,'%s%s%s','\newcommand{\WMOnum}{',float_names{I},'}'); fclose(fid);
  fid=fopen('floatsourcedef.tex','w'); fprintf(fid,'%s%s%s','\newcommand{\floatsource}{',float_dir,'}'); fclose(fid);
  fid=fopen('floatplotsdef.tex','w'); fprintf(fid,'%s%s%s','\newcommand{\floatplots}{',plot_dir,'}'); fclose(fid);
  fid=fopen('downloaddirdef.tex','w'); fprintf(fid,'%s%s%s','\newcommand{\downloaddir}{',download_dir{I},'}'); fclose(fid);


  % ------ INGEST FLOAT FILES: -----------------------------------------

  % Init of some objects:
  comments='';				% String for report on what is done here (obsolete)
  nowchar=replace(datestr(now,30),'T',''); % Datetime-stamp
  PARS=ncread(infiles{I},'PARAMETER');  % Parameters in the _prof file
  npars=size(PARS,2);			% Number of parameters in the _prof file
  % Find the complete height of matrix (not ideal method, but the
  % R-files have different length profiles and we need a matrix here):
  PRES=ncread(infiles{I},'PRES');
  m=size(PRES,1)+1; % Add one to make sure every profile ends with a NaN.
  % Find the profile-files to check here:
  Rfiles=edir(rootdirin{I});  % Will contain both D and R files, just
  %% EDIR only works on UNIX/Linux since it applies find and grep,
  % using dir here:
  %Rfiles=dir(rootdirin{I});  % Will contain both D and R files, just
                              % like on the Coriolis server.
			      % [] And A files?
  %Rfiles=strcat(rootdirin{I},{Rfiles.name}');				% Add the path
  %Rfiles=Rfiles(~contains(Rfiles,{'\.','/.'}));	% Remove the relative directory references
  % Now edir works also without unix/Linux.

  % Select file names for chosen direction:
  switch direction
   case 'A', Rfiles=Rfiles(~contains(Rfiles,'D.nc'));
   case 'D', Rfiles=Rfiles( contains(Rfiles,'D.nc')); 
  end
  n=length(Rfiles);
  % Sort the file names by cycle number in file name:
  if n>1, split(Rfiles,'_'); 
    [~,AI]=sort(ans(:,2)); Rfiles=Rfiles(AI); end
  
  if direction=='D'
    if n<1
      ['There are no direction ''',direction,''' files for float ',float_names{I},'. '];
      disp(['prepare_floats:', ans]); snippet(ans,'descending_profiles');
      continue 
    else
      ['Direction ''',direction,''' profiles have been checked for float ',...
       float_names{I},', but results are not shown in this report.'];
      snippet(ans,'descending_profiles');
    end
  end
  [PRES,PSAL,TEMP]=deal(nan(m,n)); % Originals, filled but not to be changed at all!
  [PRESqco,PSALqco,TEMPqco]=deal(repmat(' ',m,n));
  [CYCLE_NUMBER,LONG,LAT,JULD]=deal(nan(1,n));
  [POSqco,JULDqco,DIRECTION]=deal(repmat(' ',1,n));
  %[PROFILE_PSAL_QC]=deal(repmat(' ',1,n));
  [PARAMETER]=deal(repmat(' ',16,npars,1,n));
  for i=1:n							% loop Rfiles
    ncread(Rfiles{i},'LONGITUDE')';	LONG(1,i)=ans(1,1);
    ncread(Rfiles{i},'LATITUDE')';	LAT(1,i) =ans(1,1);
    ncread(Rfiles{i},'JULD')'; 		JULD(1,i)=ans(1,1);
    ncread(Rfiles{i},'PRES');		PRES(1:size(ans,1),i)=ans(:,1);
    ncread(Rfiles{i},'PSAL');		PSAL(1:size(ans,1),i)=ans(:,1);
    ncread(Rfiles{i},'TEMP');		TEMP(1:size(ans,1),i)=ans(:,1);
    ncread(Rfiles{i},'PRES_QC');	PRESqco(1:size(ans,1),i)=ans(:,1);
    ncread(Rfiles{i},'TEMP_QC');	TEMPqco(1:size(ans,1),i)=ans(:,1);
    ncread(Rfiles{i},'PSAL_QC');	PSALqco(1:size(ans,1),i)=ans(:,1);
    % ncread(Rfiles{i},'PROFILE_PRES_QC');	PROFILE_PRES_QC(1,i)=ans(1,1);
    % ncread(Rfiles{i},'PROFILE_TEMP_QC');	PROFILE_TEMP_QC(1,i)=ans(1,1);
    % ncread(Rfiles{i},'PROFILE_PSAL_QC');	PROFILE_PSAL_QC(1,i)=ans(1,1);
    ncread(Rfiles{i},'CYCLE_NUMBER')';	CYCLE_NUMBER(1,i)=ans(1,1);
    ncread(Rfiles{i},'DIRECTION')';	DIRECTION(1,i)=ans(1,1);
    ncread(Rfiles{i},'POSITION_QC')';	POSqco(1,i)=ans(1,1);
    ncread(Rfiles{i},'JULD_QC')';	JULDqco(1,i)=ans(1,1);
    ncread(Rfiles{i},'PARAMETER');
    N=size(ans,2);			% Number of parameters
    IC=ans(:,:,1,1);			% Parameter names (of first layer)
    para={'PRES','TEMP','PSAL'};	% Valid parameters for DMQC in the right order
    ic={'NB_SAMPLE_CTD','MTIME','TEMP_CNDC','PRES_MED','TEMP_MED','PSAL_MED','TEMP_STD','PSAL_STD'}; % Intermediate parameters
    for j=1:N
      pname=strip(IC(:,j)');		% The j-th parameter name
      if any(strcmp(pname,para))	% Is it a valid parameter for DMQC?
	find(strcmp(pname,para));	% Which one is it?
	PARAMETER(:,ans,1,i)=IC(:,j);	% Place it in correct column
      elseif any(strcmp(pname,ic))	% Or is it just an intermediate parameter?
	% Then deal the default calib information:
	[scientific_calib.equation.(pname)(i,:),...
	 scientific_calib.coefficient.(pname)(i,:),...
	 scientific_calib.comment.(pname)(i,:)] = deal('not applicable');
	scientific_calib.date.(pname)(i,:) = nowchar; 
	% Nothing more is supposed to be done to these intermediate parameters.
      else					% Unknown parameter!
	error([Rfiles{i},' has hitherto unknown parameter ',pname','!',...
	       ' This parameter will not get SCENTIFIC_CALIB information.',...
	       ' Check documentation and recode float ingestion in PREPARE_FLOATS.']);
      end
      % PARS is the array from the _prof file.  PARAMETER is filled from
      % the R/D-files, and its size is based on PARS.  None of these are
      % used after this, not here and not in WRITE_D.
    end % Parameter name loop
  end % R-file loop
  
  % Original number of measurements in each profile (for PROFILE_<PARAMETER>_QC in write_D):
  PROFILE_PRES_N=sum(~isnan(PRES)); PROFILE_TEMP_N=sum(~isnan(TEMP)); PROFILE_PSAL_N=sum(~isnan(PSAL));
  % Just a test on how percentage is calculated:
  % PROFILE_PSAL_QC=repmat(char(32),1,size(SALqco,2));
  % sum(SALqco=='1' | SALqco=='2' | SALqco=='5' | SALqco=='8')*100./PROFILE_PSAL_N*100;
  % PROFILE_PSAL_QC(ans==0)='F', PROFILE_PSAL_QC(ans>0)='E', PROFILE_PSAL_QC(ans>=25)='D'
  % PROFILE_PSAL_QC(ans>=50)='C', PROFILE_PSAL_QC(ans>=75)='B', PROFILE_PSAL_QC(ans==100)='A'
  
  % Get start date for datenum
  REFERENCE_DATE_TIME=ncread(Rfiles{1},'REFERENCE_DATE_TIME'); 
  % Only reading first column from file! Other profiles are not part of DMQC.

  % This part is obsolete, since descending profile(s) are removed by
  % filename above, but check in case '*D.nc' is not sufficient:
  if any(DIRECTION~=direction)
    error(['There are profiles not in the ',direction,'-direction!']);
  end
  
  % whos CYCLE_NUMBER LONG LAT JULD PRES SAL TEMP DIRECTION REFERENCE* HISTORY* 
  
  % ----- INPUT PREVIOUS QC FLAGS: -------------------------------------
  
  % Replace qco and qcn for previously DMQC'ed profiles with saved matrices: ----
  % New full length blank flag matrices to add previous DMQC to:
  [POSqcnb,JULDqcnb] = deal(repmat(' ',1,n));		
  [PRESqcnb,PSALqcnb,TEMPqcnb] = deal(repmat(' ',m,n));
  % Use qco as read from the nc-files as matrices to add saved RTQC to:
  POSqcob=POSqco; JULDqcob=JULDqco; PRESqcob=PRESqco; PSALqcob=PSALqco; TEMPqcob=TEMPqco;
  newflagfile=['savedflags_',float_names{I},'.mat'];
  if exist(newflagfile,'file')				% There has been previous DMQC.
    load(newflagfile,'*qcn','*qco');			% Load from previous DMQC.
    disp(['prepare_floats: Loading previous RTQC and DMQC flags from ',pwd,filesep,newflagfile]);
    [M,N]=size(PRESqcn);				% Size of previous flag matrices
    % Test that the saved sets of flags concur with the nc-files from
    % previous DMQC (D-files), i.e, *qco with addition of '1' and '4'
    % from *qcn should be the same as *qcob from the ncfiles:
    for i=1:length(para)					% Loop parameters
      eval(['qco=',para{i},'qco;']);				% Temporary variable for <PARAM>qco
      eval(['qco(',para{i},'qcn==''1'')=''1'';']);		% Add reversing DMQC flags to RTQC flags 
      eval(['qco(',para{i},'qcn==''4'')=''4'';']);		% Add new DMQC flags to RTQC flags 
      eval(['any(qco~=',para{i},'qcob(1:M,1:N));']);		% Compare
      if any(ans) & any(contains(Rfiles(find(ans)),'//D'))	% Warn about D-files only
	warning(['Mismatch between previously saved DMQC flags and flags in D-files for ',para{i},'!']); 
      end
    clear qco
    end   
    % Insert previous DMQC flags into the blank matrices:
    POSqcnb(1,1:N)=POSqcn; JULDqcnb(1,1:N)=JULDqcn; 	
    PRESqcnb(1:M,1:N)=PRESqcn; PSALqcnb(1:M,1:N)=PSALqcn; TEMPqcnb(1:M,1:N)=TEMPqcn;
    % Insert previous RTQC flags into the qco from nc-files matrices:
    POSqcob(1,1:N)=POSqco; JULDqcob(1,1:N)=JULDqco;
    PRESqcob(1:M,1:N)=PRESqco; PSALqcob(1:M,1:N)=PSALqco; TEMPqcob(1:M,1:N)=TEMPqco;
  end
  % Make the now filled matrices the current qcn and qco matrices:
  POSqcn=POSqcnb; JULDqcn=JULDqcnb; PRESqcn=PRESqcnb; PSALqcn=PSALqcnb; TEMPqcn=TEMPqcnb; clear *qcnb
  POSqco=POSqcob; JULDqco=JULDqcob; PRESqco=PRESqcob; PSALqco=PSALqcob; TEMPqco=TEMPqcob; clear *qcob
  % The first time around, qcn (DMQC flags) are all empty and qco (RTQC
  % flags) is as loaded from float-files. For subsequent times qcn is
  % what was done to previous profiles and blank for new profiles, and
  % qco is what has been saved for previous profiles and what has been
  % read from new (R) files. The qco matrices never change in the
  % scripts, so it will always contain RTQC only.
  %
  % Finally make temporary RTQC-flag matrices for the automatic tests
  % below to fill:
  POSqct=POSqco; JULDqct=JULDqco;
  PRESqct=PRESqco; PSALqct=PSALqco; TEMPqct=TEMPqco;
  
  % whos *qco *qct *qcn 
  
  % ---- More info for the report (including the whole appendix) only for ascending profiles: ----------------
  if direction=='A'
    % The Altimetry section
    fid=fopen('altimetry_sec.tex','w');
    if exist([download_dir{I},float_names{I},'.png'],'file')
      fprintf(fid,'%s\n','Figure~\ref{Altim} shows the comparison with altimetry.');
      fprintf(fid,'%s\n','\begin{figure}[H]');
      fprintf(fid,'%s\n','  \centering');
      fprintf(fid,'%s\n','  \includegraphics[width=\textwidth,natwidth=810,natheight=650]{\downloaddir\WMOnum .png}');
      fprintf(fid,'%s\n','  \caption{Float \WMOnum. The comparison between the sea level anomaly (SLA) from the satellite altimeter and dynamic height anomaly (DHA) extracted from the Argo float temperature and salinity. The figure is created by the CLS/Coriolis, distributed by Ifremer (ftp://ftp.ifremer.fr/ifremer/argo/etc/argo-ast9-item13-AltimeterComparison/figures/).}');
      fprintf(fid,'%s\n','  \label{Altim}');
      fprintf(fid,'%s\n','\end{figure}');
    else
      fprintf(fid,'%s\n','An altimetry report was not available at the time of reporting.');
    end
    
    % The reference data plots for the appendix:  
    ~(isnan(LONG)|isnan(LAT));wmosquares=unique(findwmo(LONG(ans),LAT(ans))); snippet(wmosquares); 
    fid=fopen('appendix.tex','w');
    for i=1:length(wmosquares)
      fprintf(fid,'%s\n','\begin{figure}[ph]');
      fprintf(fid,'%s\n','\centering');
      reffig=false;
      for j=1:length(tardir) % one fig per type of refdata
	% [] Note that page is not formatted with room for all three types. When
        % bottle data is added, will need to accomodate for this.
	wmofig=[tardir{j},filesep,tartyp{j},'_',int2str(wmosquares(i)),'.eps'];
	if exist(wmofig,'file')	% if file exist make graphics input line
	  reffig=true;
	  fprintf(fid,'%s%s%s\n','\includegraphics[width=\textwidth]{',wmofig,'}');
	end
      end 
      if reffig
	fprintf(fid,'%s%u%s\n','\caption{Overview of reference data in WMO square ',wmosquares(i),' which is traversed by Float~\WMOnum . Upper set of graphs are for CTD reference data, and lower set is for historical ARGO data (if only one set is displayed, see title of second panel). Colouring of positions in map pairs illustrate temporal coverage and latitude, respectively. The following TS and profile plots use the latter colormap.}');
      else
	fprintf(fid,'%s%u%s\n','\caption{There are no reference data in WMO square ',wmosquares(i),' traversed by Float~\WMOnum .}');
      end
      fprintf(fid,'%s%u%s\n','\label{reference',wmosquares(i),'}');
      fprintf(fid,'%s\n','\end{figure}');
    end % loop wmosquares
    fclose(fid);
  end % if direction 'A'

  % Create empty files for the report if they do not exist, as some procedures may not be performed:
  if ~exist('shipctd.tex','file'), fid=fopen('shipctd.tex','w'); fclose(fid); end
  if ~exist('coincidence.tex','file'), fid=fopen('coincidence.tex','w'); fclose(fid); end
  if ~exist('coincidence_appendix.tex','file'), fid=fopen('coincidence_appendix.tex','w'); fclose(fid); end
  if ~exist('owc.tex','file'), fid=fopen('owc.tex','w'); fclose(fid); end
  if ~exist('scientific_calib_tabular.tex','file'), fid=fopen('scientific_calib_tabular.tex','w'); fclose(fid); end
  

  % ----- Initialise the new calibration information objects: ---------- 
  % (For all profiles, but for now just with no size)
  % When data are good and no adjustment is needed and assuming that scientific_calib should be 'none' throughout:
  [scientific_calib.equation.PRES, scientific_calib.coefficient.PRES, scientific_calib.comment.PRES, ...
   scientific_calib.equation.PSAL, scientific_calib.coefficient.PSAL, scientific_calib.comment.PSAL, ...
   scientific_calib.equation.TEMP, scientific_calib.coefficient.TEMP, scientific_calib.comment.TEMP] = deal(repmat('none',n,1)); % First new entry
  % "Regardless of whether an adjustment has been applied or not, the date of delayed-mode qc for each measurement parameter should be
  % recorded in SCIENTIFIC_CALIB_DATE, in the format YYYYMMDDHHMISS."
  [scientific_calib.date.PRES		, ... 
   scientific_calib.date.PSAL		, ... 
   scientific_calib.date.TEMP	  ]	= deal(repmat(nowchar,n,1)); 
  % Known default comments when no adjustments is done:
  scientific_calib.comment.TEMP		= repmat('The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.',n,1); % Default entry
  % Default equations (PSAL is treated in WRITE_D):
  scientific_calib.equation.PRES	= repmat('PRES_ADJUSTED = PRES. ',n,1); % Default entry
  scientific_calib.equation.TEMP	= repmat('TEMP_ADJUSTED = TEMP. ',n,1); % Default entry
  % All superfluous 'none' will be removed in WRITE_D.  

  
  
  % --------- DMQC: ----------------------------------------------------------------

  
  %%%%%%% PHASE 1: Delayed-mode procedures for coordinates _and_ pressure adjustments %%%%%%%%
  %%%%%%% PHASE 2: Sea Surface Pressure adjustments %%%%%%%%
  %
  % W21 3.2. Delayed-mode procedures for JULD, LATITUDE, LONGITUDE 
  %
  % Delayed-mode operators should check that JULD in the profiles
  % are in chronological order. Erroneous or missing JULD values should
  % be replaced with another telemetered value if available, or replaced
  % with interpolated values and marked with JULD_QC = ‘8’.
  fid=fopen('JULDtest.tex','w'); [automsg,manualmsg]=deal('');
  union( find(diff(JULD)<0)+1, find(isnan(JULD)) );
  if any(ans)  
    JULD(ans)=NaN; JULDqcn(ans)='8'; JULD=efill(JULD);	% Interpolate to fill NaNs
    automsg=['Non-chronological or missing JULD in cycles ',zipnumstr(CYCLE_NUMBER(find(ans))),'.'];
    disp(['prepare_floats: ',automsg]);
    %error('Some JULD not chronological! Check!');
  end
  fprintf(fid,'%s %s\n',automsg,manualmsg);
  fclose(fid);
  % Calculate DATES as OWC wants them:
  REFERENCE_DATE_TIME'; time=datenum([ans(1:8),'T',ans(9:14)],'yyyymmddTHHMMSS')+JULD; 
  DATES=dyear(time);
  [~,MONTHS]=datevec(time); % For the seasonal aspect
  % JULDqcn is set now.
  
  % Sort cycle number (not sure this should be done, so just a check and error for now):
  % Double cycle numbers messes up the rest. Pick out upcasts.
  % Also Rfiles may be ordered wrong.
  [ans,IA]=sort(CYCLE_NUMBER); % C = A(IA);
  if ~all(ans==CYCLE_NUMBER)
    error('Cycle numbers not in temporal succession!');
    LONG=LONG(IA); LAT=LAT(IA); DATES=DATES(IA); CYCLE_NUMBER=CYCLE_NUMBER(IA); PRES=PRES(:,IA); PSAL=PSAL(:,IA); TEMP=TEMP(:,IA);
    CYCLE_NUMBER=ans;
  end
  snippet([int2str(CYCLE_NUMBER(1)),'-',int2str(CYCLE_NUMBER(end))],'cycle-numbers');
  
  % Size of data geometry is final after sorting:
  [m,n]=size(PRES);
  PROFILE_NO=1:n;
  % Default adjusted values (after sorting):
  PRES_ADJUSTED=PRES; TEMP_ADJUSTED=TEMP; PSAL_ADJUSTED=PSAL;

  % [] How to know about wrong dates not being just shifted? Relation
  % to cycle number?

  % Profile positions in LONGITUDE, LATITUDE should be checked for
  % outliers.  Erroneous or missing LONGITUDE, LATITUDE values should
  % be replaced with another telemetered value if available, or replaced
  % with interpolated values and marked with POSITION_QC = ‘8’.
  %
  % Automated pre-check: 
  if any(PROFILE_NO>checked{I}) 
    fid=fopen('POStest.tex','w'); [automsg,manualmsg]=deal('');
    % Extrapolated data, i.e. 'interpolated' data at the end of series
    % cannot be used (extrapolation makes no sense):
    j=find(POSqco=='8' | POSqco=='9');
    if any(j), if j(end)==n			% Missing position at end 
      [~,gi]=groups(diff(j)); 
      if j(gi(end))+1==n		% Several, actually
	j=j(gi(1,end):gi(2,end)+1);
	POSqct(j)='4';			% Flag them
	automsg=['There were extrapolated (!) or missing positions at the end of series. ' ...
		 'Hence, positions of cycles ',zipnumstr(CYCLE_NUMBER(j)),' as a rule flagged as bad.'];
	disp(['prepare_floats: ',automsg]);
      end
    end, end
    LAT<-90 | 90<LAT;
    if any(ans)
      POSqct(ans)=='4';
      'There were latitudes outside 90S and 90N.'; disp(['prepare_floats: ',ans]); automsg=[automsg,' ',ans];
    end
    % Visual test of positions:
    POSqct = check_profiles(CYCLE_NUMBER,LONG,LAT,POSqct);
    % Add only different flags as new flags:
    POSqct~=POSqco; POSqcn(ans)=POSqct(ans);	
    if any(ans)
      manualmsg=['Automatic and visual DMQC have changed the position flags of cycles ',zipnumstr(CYCLE_NUMBER(find(ans))),'.'];
      disp(['prepare_floats: ',manualmsg]);
    end
    fprintf(fid,'%s %s\n',automsg,manualmsg);
    fclose(fid);
  end
  
  % Apply LATITUDE, LONGITUDE flags immediately in order for the visual tests and plots hereafter to function: 
  POSqco=='4' & POSqcn~='1' | POSqcn=='4'; LONG(ans)=NaN; LAT(ans)=NaN; 
  
  % whos LONG LAT POSqcn JULDqcn 
  % whos *_ADJUSTED
  
  
    if direction=='D'
      % Apply salinity flags from the RTQC in the surface, to the
      % _ADJUSTED already here since the descending profiles contain
      % above surface salinities. This is done mainly to avoid too many
      % unneccessary detections in the automatic test below, and not to
      % zoom out the automatic profile plots with far off values. This
      % does not affect the visual control, as that uses non adjusted
      % data.
      PSALqco~='1' & PRES<10;		% Find data to omit
      igi=find(ans);			% For testing after
      igiS=PSAL_ADJUSTED(igi);		% Save to put back in after the tests and checks
      PSAL_ADJUSTED(ans)=NaN;		% Omit data for the tests
      PSALqct(ans & PSALqco~='4')='Y';	% Assign new flags (used in manual control part)
      PSALqcn(ans & PSALqco~='4')='4';	% Assign new flags (used when no manual control)
      % Special case: Just omit and flag all shallow (useless) data:
      % if contains(platyp,{'PROVOR'})
      % 	PRES<10;
      % 	TEMPqct(ans)='Y';	% Assign new flags (used in manual control part)
      % 	PSALqct(ans)='Y';	% Assign new flags (used in manual control part)
      % end 
      %
      % Skipped, because it looks not so good to do this
      % indiscriminately. The user must just realise that descending
      % floats are not trustworthy in the top.
    end

    
    
    
  % Plot the input data raw:
  DENS = sw_pden(PSAL_ADJUSTED,TEMP_ADJUSTED,PRES_ADJUSTED,0);	% Potential density
  zerow=zeros(1,n);						% Row of zeros
  plot_profiles;						% Plots based on ADJUSTED parameters
  print(gcf,'-depsc',[outfiles{I}(1:end-4),'_raw.eps']);	% Print to file
  
  % whos DENS

  % W21 3.3. Delayed-mode procedures for pressure
  %
  % Bad data points identified by visual inspection from delayed-mode
  % analysts are recorded with PRES_ADJUSTED_QC = ‘4’ and PRES_QC =
  % '4'. For these bad data points, TEMP_QC, TEMP_ADJUSTED_QC, PSAL_QC,
  % PSAL_ADJUSTED_QC should also be set to ‘4’.  Please note that
  % whenever PARAM_ADJUSTED_QC = ‘4’, both PARAM_ADJUSTED and
  % PARAM_ADJUSTED_ERROR should be set to FillValue.
  %
  % √ This is taken care of both here for OWC and in WRITE_D.
  %
  % PRESqcn TEMPqcn PSALqcn % Change QC flags if necessary, all in this case.
  %
  % Check for negative pressures:
  PRES_ADJUSTED<0;
  if any(ans,'all'), 
    warning('Negative pressures detected prior to correction!'); 
  end
  % Pressure adjustments for APEX, or not:
  fid=fopen('pressure-adjustment.tex','w');
  fprintf(fid,'%s\n','Sea surface pressure adjustments should be done for APEX floats (Wong et al., 2021). ');
  if contains(platyp,'APEX') 
    % 3.3.1 Delayed-mode pressure adjustment for APEX and NAVIS floats
    SP.CN=ncread(techfiles{I},'CYCLE_NUMBER')';
    SP.TPN=string(ncread(techfiles{I},'TECHNICAL_PARAMETER_NAME')');
    SP.TPV=ncread(techfiles{I},'TECHNICAL_PARAMETER_VALUE')';
    iSP=find(contains(SP.TPN,'PRES_SurfaceOffset'));
    if ~isempty(iSP)
      fprintf(fid,'%s%s%s%s\n','This is an ',...
	      platyp,...
	      ' float. Upper panel of Figure~\ref{fig:pres} shows the original and despiked PRES$\_$SurfaceOffsetNotTruncated$\_$dbar (SP)', ...
	      ' in addition to the adjusted pressures from the top of each profile. Adjustment is done when valid SP values are available. ');
      SP.CN=SP.CN(iSP); SP.TPN=SP.TPN(iSP); SP.TPV=str2num(SP.TPV(iSP,:));
      % (1):
      unique(SP.TPN);
      if length(ans)>1, error('Technical parameter name for SP not unique!'); end
      if contains(ans,'PRES_SurfaceOffsetTruncatedPlus5dbar_dbar')
	SP.TPV=SP.TPV-5; 
	error('Truncated surface pressure! You need to code for TNPD.');
      end
      % (2):
      diff(SP.TPV); SP.TPV(find((ans(1:end-1)>5 & ans(2:end)<-5) | (ans(1:end-1)<-5 & ans(2:end)>5))+1)=NaN;
      SP.mn=mean(SP.TPV,'omitnan');
      SP.mf=medfilt1(SP.TPV-SP.mn,5)+SP.mn;
      SP.SP=SP.TPV;
      SP.SP(abs(SP.SP-SP.mf)>1)=NaN;
      % (3):
      SP.SP=efill(SP.SP)';
      % (4):
      % Get Cond to re-calculate PSAL:
      Cond=nan(m,n); 
      ii=find(~isnan(PRES)&~isnan(PSAL)&~isnan(TEMP));		% All parameters non-NaN because 
      Cond(ii)=gsw_C_from_SP(PSAL(ii),TEMP(ii),PRES(ii));	% gsw_ can't take any NaNs! :-(
      % PRES_ADJUSTED (cycle i) = PRES (cycle i) – SP (cycle i+1).
      [C,SP.IA,SP.IB] = intersect(CYCLE_NUMBER+1,SP.CN,'stable'); % IA and IB such that C = A(IA,:) and C = B(IB,:).
      PRES_ADJUSTED(:,SP.IA) = PRES_ADJUSTED(:,SP.IA) - SP.SP(SP.IB)';
      PSAL_ADJUSTED(ii) = gsw_SP_from_C(Cond(ii),TEMP_ADJUSTED(ii),PRES_ADJUSTED(ii)); % SP = gsw_SP_from_C(C,t,p)	% Re-calculate SAL

      % Calibration info:
      'PRES_ADJUSTED = PRES - dP. '; 
		scientific_calib.equation.PRES(SP.IA,1:size(ans,2)) = repmat(ans,length(SP.IA),1);  % First DMQC entry; overwrite
      char(strcat({'dP = '},num2str(SP.SP(SP.IB)),' dbar. '));
		scientific_calib.coefficient.PRES(SP.IA,1:size(ans,2)) = ans;  % First DMQC entry; overwrite
      'Pressures adjusted by using pressure offset at the sea surface. The quoted error is manufacturer specified accuracy in dbar. ';
		scientific_calib.comment.PRES(SP.IA,1:size(ans,2)) = repmat(ans,length(SP.IA),1); % First DMQC entry; overwrite
		scientific_calib.date.PRES(SP.IA,:)	= repmat(nowchar,length(SP.IA),1); % Always overwrite date      

      % 3.3.2. Truncated negative pressure drift (TNPD) in APEX floats
      % [] Not relevant for our NorARGO floats, so not implemented. Sorry! Feel free to add code and notify on github.
    else
      clear SP
      fprintf(fid,'%s%s%s%s\n','This is an ',...
	      platyp,...
	      ' float, but it has no valid PRES$\_$SurfaceOffsetNotTruncated$\_$dbar (SP) values. ',...
	      ' Instead, upper panel of Figure~\ref{fig:pres} shows the surface pressure from the top of each profile. ');
    end
  else % Ikke APEX-float, men plot tid mot overflatetrykk, sensortrykk.
    fprintf(fid,'%s%s%s\n','This is an ',...
	    platyp,...
	    ' float. Instead, upper panel of Figure~\ref{fig:pres} shows the surface pressure from the top of each profile. ');
  end
  fprintf(fid,'%s\n','\begin{figure}[H]');
  fprintf(fid,'%s\n','\centering');
  fprintf(fid,'%s\n','\includegraphics[width=\textwidth,natwidth=1500,natheight=1250]{\floatsource\WMOnum_PRES.png}');
  fprintf(fid,'%s\n',['\caption{Float \WMOnum\ pressure data. ',...
		      'Upper panel: Top of profile pressure series. ',...
		      'Lower panel: Pressure data for all profiles presented by every 10th (and the deepest) measurement in each profile. ',...
		      'Blue dots indicate pressure value in the real-time. ',...
		      'Any extra marks show RTQC and DMQC flags as explained in the legend. ',...
		      'Any ''4'', or ''8'' are shown regardless of depth, not just at every 10th.}']);
  fprintf(fid,'%s\n','\label{fig:pres}');
  fprintf(fid,'%s\n','\end{figure}');
  fclose(fid);
  % 'PRES_ADJUSTED = PRES - dP'; 

  % Check for negative pressures after the correction:
  PRES_ADJUSTED<0;
  if any(ans,'all')
    warning('Negative pressures detected after the correction! (removed)');
    PRESqct(ans)='X'; % For the visual control below
  end  
  
  % whos PRES_ADJUSTED PSAL_ADJUSTED PRESqct
    
  % For the plots later, but also PTMPo here:
  PTMP = sw_ptmp(PSAL,TEMP,PRES,0); 
  
  % whos PTMP      

  % Save ASD profiles salinities and ptmp from PSAL and PTMP (before test plots):
  ASD=n+1;			% Add 1 just since -1 is used later
  % For stability, extract your own explicitly set calseries in set_calseries.m.
  jj=find(cal_action{I}==4);					% Which cal action (if any) is 4? 
  aa=strip(readlines([my_working_dir,'DMQC',filesep,float_names{I},filesep,'set_calseries.m']));	% Read local set_calseries.m
  find(contains(aa,'calseries = [ones') & ~startsWith(aa,'%')); % Scan set_calseries.m 
  eval(aa(ans(end)));						% Get the last definition of calseries
  if any(jj) & direction=='A' & exist('calseries','var')
    [C,IA] = unique(calseries,'stable');		% Find the shifts in calseries
    ASD=IA(jj);						% The one before the first profile of the jjth action (4). 
    SALo=PSAL(:,ASD:end);				% Save to fill SAL with in the end
    PTMPo=PTMP(:,ASD:end);				% Save to fill PTMP with in the end
    disp(['prepare_floats: Conserving salinity profiles with ASD (profiles ',int2str(ASD),'-',int2str(n),').']);
  end

  % whos ASD SALo PTMPo
  
  
  % ----- PHASE 3 and PHASE 4 ------------------------------------------
  
  % W21 3.4. Delayed-mode procedures for temperature
  %
  % Bad data points identified by visual inspection from delayed-mode
  % analysts are recorded with TEMP_ADJUSTED_QC = ‘4’ and TEMP_QC =
  % '4'. Please note that whenever PARAM_ADJUSTED_QC = ‘4’,
  % PARAM_ADJUSTED = FillValue, and PARAM_ADJUSTED_ERROR = FillValue.
  %
  % TEMP_ADJUSTED, TEMP_ADJUSTED_ERROR, and TEMP_ADJUSTED_QC should be
  % filled even when the data are good and no adjustment is needed. In
  % these cases, TEMP_ADJUSTED_ERROR can be the manufacturer’s quoted
  % accuracy at deployment, which is 0.002°C.
  % √ Done in WRITE_D.
  %
  % Please use the SCIENTIFIC CALIBRATION section in the netCDF files to
  % record details of the delayed-mode adjustment.
  %
  % TEMPqc % Change QC flags if necessary
	
  % W21 3.5. Delayed-mode procedures for salinity
  %
  % It is recommended that float salinity be adjusted for pressure
  % offset and cell thermal mass error before sensor drift
  % adjustment.
  % [√] Done above.
  %
  %W Operators should also ensure that other float measurements (PRES,
  %W TEMP, LATITUDE, LONGITUDE, JULD) are accurate or adjusted before
  %W they input them into the statistical tools for estimating reference
  %W salinity.
  %
  % PSALqcn % Change QC flags if necessary

  J=find(PROFILE_NO>checked{I});
  %J=[]; visd(2)=0; % Use this for no DMQC checks
  
  if any(J)	% IF NEW PROFILES DO AUTOMATIC AND MANUAL CONTROL OF _NEW_PROFILES_
  
    %%%%%%% PHASE 3: Automated tests prior to visual control %%%%%%%

    % These tests:
    %	- Checks new/wanted profiles
    %	- Plots warning plots not used in report
    %	- Flags with different letters in temprary flag matrices
    
    % NOTE: For Deep Arvor profiles, values deeper than 2000 m PRES and
    % TEMP are flagged '2' and PSAL flagged '3' already in R-files. The
    % same is the case for previously DMQC'd files (D-files).
 
    % Initial check, since the code here adds flag 'X' to errors found in automatic tests: 
    if any(PRESqco=='Y' | PSALqco=='Y' | TEMPqco=='Y' | PRESqcn=='Y' | PSALqcn=='Y' | TEMPqcn=='Y' | ...
	   PRESqco=='X' | PSALqco=='X' | TEMPqco=='X' | PRESqcn=='X' | PSALqcn=='X' | TEMPqcn=='X' | ...
	   PRESqco=='S' | PSALqco=='S' | TEMPqco=='S' | PRESqcn=='S' | PSALqcn=='S' | TEMPqcn=='S' | ...
	   PRESqco=='D' | PSALqco=='D' | TEMPqco=='D' | PRESqcn=='D' | PSALqcn=='D' | TEMPqcn=='D' | ...
	   PRESqco=='G' | PSALqco=='G' | TEMPqco=='G' | PRESqcn=='G' | PSALqcn=='G' | TEMPqcn=='G' | ...
	   PRESqco=='I' | PSALqco=='I' | TEMPqco=='I' | PRESqcn=='I' | PSALqcn=='I' | TEMPqcn=='I')
      error('There are already flag Y, X, S, D, G or I from RTQC or previous DMQC! Must recode PREPARE_FLOATS.'); 
    end
    
    DENS = sw_pden(PSAL_ADJUSTED,TEMP_ADJUSTED,PRES_ADJUSTED,0);	% For the plotting of warning plots
    zerow=zeros(1,n);							% Row of zeros used in tests
    [snb2,tnb2,snb1,tnb1,snbG,tnbG]=deal([]);
    
    % Pressure increasing test / monotonically increasing pressure test:
    logical([zerow;diff(PRES_ADJUSTED,1,1)<=0]);  % The testvalue
    ans(PRES_ADJUSTED<=MAP_P_EXCLUDE)=logical(0); % Do not be concerned with monotonicity above this depth.
    if any(ans,'all')
      jnb=find(any(ans)); pnb=ans; nbt='Non-monotonic pressure';
      disp(['prepare_floats: ',nbt,'in float ',float_names{I},' cycle(s): ',int2str(CYCLE_NUMBER(jnb)),' !']);	% Only warning
      plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_pressure_warning.eps'])
      % [] If there is a region of constant pressure, all but the
      % first of the consecutive levels of constant pressure should be
      % flagged as bad data (‘4’).
      % If there is a region where pressure reverses, all of the
      % pressures in the reversed part of the profile should be flagged
      % as bad data. 
      % All pressures flagged as bad data and all of the associated
      % temperatures and salinities should be removed.
      % NO ACTION. IMPLEMENT WHEN IT HAPPENS.
      % Monotonicity is a test in CHECK_PROFILES below.
      comments=[comments,lower(nbt),' (>',int2str(MAP_P_EXCLUDE),' m); '];
      clear jnb nbt
    else
      system(['rm -f ',outfiles{I}(1:end-4),'_pressure_warning.eps']); % Clean up any previous plot
    end
    
    % NEW! THE DOUBLE-POINTED SPIKE TEST:
    % Some times the spike is made up of two deviating points.
    % Test value = | (V2 + V3)/2 – (V1 + V4)/2 | – | (V3 – V2) / 2 | – | (V4 – V1) / 2 |
    % For this test the testvalues are placed on the top point of the
    % imagined double point spike, hence, any finds will have to be
    % supplemented with a true value below as well. There is no chance
    % of looping over to next profile, since we pad with zeros at the
    % bottom of the testvalue matrix. We use the same criteria as for
    % single spikes.
    V=PSAL_ADJUSTED;
    testvalue = [zerow; abs( (V(2:end-2,:)+V(3:end-1,:))/2   -      (V(4:end,:)+V(1:end-3,:))/2 ) ...
		      - abs( (V(2:end-2,:)-V(3:end-1,:))/2 ) - abs( (V(4:end,:)-V(1:end-3,:))/2 ) ; zerow ; zerow];
    [PRES_ADJUSTED<500 & testvalue>0.9 | PRES_ADJUSTED>=500 & testvalue>0.02 | PRES_ADJUSTED>=1000 & testvalue>0.005]; % By experience clear spikes in the Nordic Seas
    ans(find(ans)+1)=logical(1); % Also the next value
    if any(ans,'all')
      jnb=find(any(ans)); snb2=ans; nbt='Double-pointed salinity spike(s)';
      disp(['prepare_floats: ',nbt,' in float ',float_names{I},' cycle(s): ',int2str(CYCLE_NUMBER(jnb)),' !']);
      plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_S-double_spike_warning.eps']);
      PSALqct(snb2)='D'; 
      comments=[comments,lower(nbt),'; '];
      clear jnb nbt
    else
      system(['rm -f ',outfiles{I}(1:end-4),'_S-double_spike_warning.eps']);
    end
    V=TEMP_ADJUSTED;
    testvalue = [zerow; abs( (V(2:end-2,:)+V(3:end-1,:))/2   -      (V(4:end,:)+V(1:end-3,:))/2 ) ...
		      - abs( (V(2:end-2,:)-V(3:end-1,:))/2 ) - abs( (V(4:end,:)-V(1:end-3,:))/2 ) ; zerow ; zerow];
    [PRES_ADJUSTED<500 & testvalue>6 | PRES_ADJUSTED>=500 & testvalue>2 | PRES_ADJUSTED>=1000 & testvalue>1]; % By experience in the Nordic Seas  
    ans(find(ans)+1)=logical(1); % Also the next value
    if any(ans,'all')
      jnb=find(any(ans)); tnb2=ans; nbt='Double-pointed temperature spike(s)';
      disp(['prepare_floats: ',nbt,' in float ',float_names{I},' cycle(s): ',int2str(CYCLE_NUMBER(jnb)),' !']);
      plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_T-double_spike_warning.eps']);
      TEMPqct(tnb2)='D';
      comments=[comments,lower(nbt),'; '];
      clear jnb nbt
    else
      system(['rm -f ',outfiles{I}(1:end-4),'_T-double_spike_warning.eps']);
    end
    
    % Spike tests (RTQC double check and stricter DMQC check for some):
    % Test value = | V2 – (V3 + V1)/2 | – | (V3 – V1) / 2 |
    % according to EuroGOOS,  where V2 is the measurement being tested
    % as a spike, and V1 and V3 are the values above and below. 
    V=PSAL_ADJUSTED; V(snb2)=NaN;	% Remove above test results
    testvalue = [zerow; abs(V(2:end-1,:)-(V(3:end,:)+V(1:end-2,:))/2) - abs((V(3:end,:)-V(1:end-2,:))/2) ; zerow];
    [PRES_ADJUSTED<500 & testvalue>0.9 | PRES_ADJUSTED>=500 & testvalue>0.02 | PRES_ADJUSTED>=1000 & testvalue>0.005]; % By experience clear spikes in the Nordic Seas
    if any(ans,'all')
      jnb=find(any(ans)); snb1=ans; nbt='Salinity spike(s)';
      disp(['prepare_floats: ',nbt,' in float ',float_names{I},' cycle(s): ',int2str(CYCLE_NUMBER(jnb)),' !']);
      plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_S-spike_warning.eps']);
      PSALqct(snb1)='S'; 
      comments=[comments,lower(nbt),'; '];
      clear jnb nbt
    else
      system(['rm -f ',outfiles{I}(1:end-4),'_S-spike_warning.eps']);
    end
    V=TEMP_ADJUSTED; V(tnb2)=NaN;	% Remove above test results
    testvalue = [zerow; abs(V(2:end-1,:)-(V(3:end,:)+V(1:end-2,:))/2) - abs((V(3:end,:)-V(1:end-2,:))/2) ; zerow];
    [PRES_ADJUSTED<500 & testvalue>6 | PRES_ADJUSTED>=500 & testvalue>2 | PRES_ADJUSTED>=1000 & testvalue>1]; % By experience in the Nordic Seas  
    if any(ans,'all')
      jnb=find(any(ans)); tnb1=ans; nbt='Temperature spike(s)';
      disp(['prepare_floats: ',nbt,' in float ',float_names{I},' cycle(s): ',int2str(CYCLE_NUMBER(jnb)),' !']);
      plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_T-spike_warning.eps']);
      TEMPqct(tnb1)='S';
      comments=[comments,lower(nbt),'; '];
      clear jnb nbt
    else
      system(['rm -f ',outfiles{I}(1:end-4),'_T-spike_warning.eps']);
    end
  
    % Gradient test:
    % Test value = | V2 − (V3 + V1)/2 |
    % where V2 is the measurement being tested, and V1 and V3 are the values above and below.
    V=PSAL_ADJUSTED; V(snb2)=NaN;  V(snb1)=NaN; % Remove above test results
    testvalue = [zerow ; abs(V(2:end-1,:)-(V(3:end,:)+V(1:end-2,:))/2) ; zerow];
    switch direction
     case 'A'
      [PRES_ADJUSTED<500 & testvalue>1.5 | PRES_ADJUSTED>=500 & testvalue>0.5]; % The RTQC9 limits 
     case 'D'
      % The RTQC9 limits plus ignore near surface which is often hampered by bad salinities in descending profiles 
      [PRES_ADJUSTED>2 & (PRES_ADJUSTED<500 & testvalue>1.5 | PRES_ADJUSTED>=500 & testvalue>0.5)]; 
    end
    if any(ans,'all')
      jnb=find(any(ans));  snbG=ans; nbt='Salinity gradient(s)';
      disp(['prepare_floats: ',nbt,' in float ',float_names{I},' cycle(s): ',int2str(CYCLE_NUMBER(jnb)),' !']);
      plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_S-gradient_warning.eps']);
      PSALqct(snbG)='G';
      comments=[comments,lower(nbt),'; '];
      clear jnb nbt
    else
      system(['rm -f ',outfiles{I}(1:end-4),'_S-gradient_warning.eps']);
    end
    V=TEMP_ADJUSTED; V(tnb2)=NaN;  V(tnb1)=NaN;	% Remove above test results
    testvalue = [zerow ; abs(V(2:end-1,:)-(V(3:end,:)+V(1:end-2,:))/2) ; zerow];
    [PRES_ADJUSTED<500 & testvalue>9 |  PRES_ADJUSTED>=500 & testvalue>3]; % The RTQC9 limits 
    if any(ans,'all')
      jnb=find(any(ans)); tnbG=ans; nbt='Temperature gradient(s)';
      disp(['prepare_floats: ',nbt,' in float ',float_names{I},' cycle(s): ',int2str(CYCLE_NUMBER(jnb)),' !']);
      plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_T-gradient_warning.eps']);
      TEMPqct(tnbG)='G';
      comments=[comments,lower(nbt),'; '];
      clear jnb nbt 
    else
      system(['rm -f ',outfiles{I}(1:end-4),'_T-gradient_warning.eps']);
    end

    % Density inversion test:
    pres=PRES_ADJUSTED; pres(PRESqco=='4'|PSALqco=='4'|TEMPqco=='4')=NaN; % Skip already flagged data
    PR = pres - [zerow; nandiff(pres)/2];		% Reference pressure
    dens = sw_pden(PSAL_ADJUSTED,TEMP_ADJUSTED,pres,PR);% Potential density
    dens(snb2)=NaN; dens(tnb2)=NaN; dens(snb1)=NaN; dens(tnb1)=NaN; dens(snbG)=NaN; dens(tnbG)=NaN; % Remove above test results
    downw=nandiff(dens)<-0.03;				% Downward test
    pres([downw;logical(zerow)])=NaN;			% Remove found
    PR = pres - [zerow; nandiff(pres)/2];		% Reference pressure again
    dens = sw_pden(PSAL_ADJUSTED,TEMP_ADJUSTED,pres,PR);% Potential density again
    upw=nandiff(dens)<-0.03;				% Upward test
    logical([downw;zerow]+[zerow;upw]);			% Logical for bad data
    if any(ans,'all')
      jnb=find(any(ans)); inb=ans; nbt='Density inversion(s)';
      disp(['prepare_floats: ',nbt,' in float ',float_names{I},' cycle(s): ',int2str(CYCLE_NUMBER(jnb)),' !']);	% Only warning
      plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_inversion_warning.eps'])
      PSALqct(inb)='I'; TEMPqct(inb)='I';
      comments=[comments,lower(nbt),'; '];
      clear inb jnb nbt
    else
      system(['rm -f ',outfiles{I}(1:end-4),'_inversion_warning.eps']);
    end

    % A little text about profiles with RTQC flags (in all profiles):
    jnb=find(any(PRESqco=='4'|TEMPqco=='4'|PSALqco=='4'));
    if any(jnb)
      disp(['prepare_floats: RTQC flags ''4'' exists in float ',float_names{I},' cycle(s): ',int2str(CYCLE_NUMBER(jnb)),' !']);
    else
      disp(['prepare_floats: No RTQC flags ''4'' exists in float ',float_names{I},'.']);
    end
    clear jnb

    
    %%%%%%% PHASE 4: Visual control of all profiles and flags %%%%%%%

    % Load the oposite direction profiles to use as parallel profiles in
    % visual control:
    switch direction
     case 'D', [outfiles{I}(1:end-5),outfiles{I}(end-3:end)];     % Remove the D
     case 'A', [outfiles{I}(1:end-4),'D',outfiles{I}(end-3:end)]; % Add the D
    end
    try
      pp=load(ans,'PRES','TEMP','PSAL','CYCLE_NUMBER');
      % In order to send these to check_profiles as the matching oposite
      % parallel profiles, the cycle numbers have to be matched up:
      jj=find(ismember(pp.CYCLE_NUMBER,CYCLE_NUMBER));
      pp.CYCLE_NUMBER=pp.CYCLE_NUMBER(jj); 
      pp.PRES=pp.PRES(:,jj); pp.TEMP=pp.TEMP(:,jj); pp.PSAL=pp.PSAL(:,jj);
    catch, pp=[];
    end
    % If the number of parallel profiles is less than the current
    % dataset, check_profiles handles it.  In the automatic tests,
    % temporary flags are given, but _ADJUSTED variables are NOT changed
    % (NaN'ed) anymore. Also, the flags given there are qct='X', so that
    % they can be visually inspected and marked below. The warnings from
    % the above tests should be duly noted while doing the visual check.
    % Apart from that, at this point qct is identical to qco (apart from
    % the top-removals in descending profiles), and qcn is flagging from
    % previous DMQC.
    %
    [PRESqct,TEMPqct,PSALqct] = check_profiles(PRES,TEMP,PSAL,PRESqct,TEMPqct,PSALqct,PROFILE_NO(J),LONG,LAT,CYCLE_NUMBER,'refdatadir',tardir,'time',time,'parallelprofiles',pp);
    %
    % Change any flags from the automatic checks overlooked in visual:
    PRESqct(PRESqct=='Y')='4'; TEMPqct(TEMPqct=='Y')='4'; PSALqct(PSALqct=='Y')='4';
    PRESqct(PRESqct=='X')='4'; TEMPqct(TEMPqct=='X')='4'; PSALqct(PSALqct=='X')='4';
    PRESqct(PRESqct=='S')='4'; TEMPqct(TEMPqct=='S')='4'; PSALqct(PSALqct=='S')='4';
    PRESqct(PRESqct=='D')='4'; TEMPqct(TEMPqct=='D')='4'; PSALqct(PSALqct=='D')='4';
    PRESqct(PRESqct=='G')='4'; TEMPqct(TEMPqct=='G')='4'; PSALqct(PSALqct=='G')='4';
    PRESqct(PRESqct=='I')='4'; TEMPqct(TEMPqct=='I')='4'; PSALqct(PSALqct=='I')='4';
    % If you started the visual not at the first profile or stopped the
    % visual before checking all profiles, the remaining columns of qct
    % will be blank. Automatic flags from those columns will be lost
    % too, since they are not controlled visually. The first and last
    % checked profile will then be number:
    find((~all(PRESqct==' ' & TEMPqct==' ' & PSALqct==' ')));
    visd=mima(ans);
    % Now change new flags only when 
    % - profile is checked in present DMQC
    % - and temporary flags differ from RTQC flags
    mask=false(m,n); mask(:,visd(1):visd(2))=true;
    PRESqcn(mask)=' '; TEMPqcn(mask)=' '; PSALqcn(mask)=' ';	% Clear qcn for checked profiles  
    PRESqct~=PRESqco & mask ; PRESqcn(ans)=PRESqct(ans);	% Add only new p flags in checked profiles
    TEMPqct~=TEMPqco & mask ; TEMPqcn(ans)=TEMPqct(ans);	% Add only new T flags in checked profiles
    PSALqct~=PSALqco & mask ; PSALqcn(ans)=PSALqct(ans);	% Add only new S flags in checked profiles
    % % - temporary flags differ from old (RTQC or previous DMQC) flags,
    % %   i.e., new flags or reversal of RTQC flags; 
    % % - and temporary flags differ from DMQC flags, i.e. DMQC is
    % %   re-done on some profiles (changed) or on new profiles (from
    % %   blanks); 
    % % - but only overwrite profiles that you actually have checked
    % %   this time. 
    % % mask=false(m,n); mask(:,visd(1):visd(2))=true; 
    % % PRESqct~=PRESqco & PRESqct~=PRESqcn & mask; PRESqcn(ans)=PRESqct(ans);
    % % TEMPqct~=TEMPqco & TEMPqct~=TEMPqcn & mask; TEMPqcn(ans)=TEMPqct(ans);
    % % PSALqct~=PSALqco & PSALqct~=PSALqcn & mask; PSALqcn(ans)=PSALqct(ans);
    %clear *qct
    % Reversal of old DMQC is also possible and fine, since qcn is
    % updated and saved after each session.
    % qct is now gone, while *qco and *qcn do not change from here on until saved.
    %
    disp(['prepare_floats: Profiles ',int2str(visd(1)),'-',int2str(visd(2)),' checked. Update the parameter ''',...
	  checked_name,''' in INIT_DMQC from ',int2str(checked{I}),' to ',int2str(visd(2)),...
	  ' for Float ',float_names{I},', and re-run prepare_floats!']);

    
  else		% IF NO NEW/UNCHECKED PROFILES => MAKE THE MULTIPANEL PLOTS: 
		% This is the case when visual control is already done on all
		% profiles, but you re-run this script in order to make the
		% multipanel plots about the flags. 
		
    fid=fopen('RTQCcheck.tex','w'); % For the report
    
    % Find profiles with any bad-flags:
    badpoints = PRESqco=='4' | PSALqco=='4' | TEMPqco=='4' | PRESqcn=='4' | PSALqcn=='4' | TEMPqcn=='4';
    J=find(any(badpoints(:,1:ASD-1),1)); % Limit to cycles prior to ASD
    % If any bad flags, make plots, but only for ascending profiles (descending almost always have flags)
    if any(J) & direction=='A'
      N=length(J);
      for j=1:N		% Loop all columns with bad-flags
	% Vertical indices for old and new flags in this profile:
	iPo=find(PRESqco(:,J(j))=='4'); iPn=find(PRESqcn(:,J(j))=='4');
	iTo=find(TEMPqco(:,J(j))=='4'); iTn=find(TEMPqcn(:,J(j))=='4');
	iSo=find(PSALqco(:,J(j))=='4'); iSn=find(PSALqcn(:,J(j))=='4');
	jj=J(j)+[-3:3];					% indices for neighbouring profiles
	jj=jj(0<jj&jj<=n);				% keep indices inside range 
	~isnan(LONG(jj))&~isnan(LAT(jj)); jj=jj(ans);	% Avoid NaNs (POSqco is not applied yet)
	
	if ~any(ismember(J(j),jj))		% profile not in valid 'neighbourhood'
	  jj=J(j);				% so just use profile position
	end
	if ismember(j,1:5:N)				% Make new figure window
	  ju=j-1;					% j used in previous figure(s)
	  M=min(5,N-ju);				% Rows in this figure
	  figure; shape=[1 1 1000 200*M];
	  set(gcf,'units','points','innerposition',shape,...
		  'paperunits','points','paperposition',shape,'PaperSize',shape(3:4),...
		  'PaperPositionMode','manual','RendererMode','manual','Renderer','opengl');
	end
	% Referencedata:
	if ~all(isnan(LONG(jj))) | ~all(isnan(LAT(jj)))
	  [LO{j},LA{j}]=halo(LONG(jj),LAT(jj),.3,.2);
	  [pres,temp,sal,long,lat,~,~,dt,dmon]=inpolygon_referencedata(LO{j},LA{j},tardir,[],time(jj));
	  dtnote={'Reference data''s seasonal';['offset is ',int2str(min(dmon)),'-',int2str(max(dmon)),' months']};
	else 
	  [pres,temp,sal,long,lat,~,dt,dmon,LO{j},LA{j}]=deal(NaN);
	  dtnote='';
	end
	% Temperature plot:
	aT=subplot(M,2,2*((j-ju)-1)+1); 
	hrT=plot(temp,pres,'color',grey); 
	hT=line(TEMP(:,jj),PRES(:,jj)); heT=line(TEMP(:,J(j)),PRES(:,J(j))); 
	if any(iTo), hEp=line(TEMP(iTo,J(j)),PRES(iTo,J(j))); set(hEp,'color','m','marker','s','linestyle','none'); end
	if any(iPo), hEp=line(TEMP(iPo,J(j)),PRES(iPo,J(j))); set(hEp,'color','g','marker','s','linestyle','none'); end
	if any(iTn), hEp=line(TEMP(iTn,J(j)),PRES(iTn,J(j))); set(hEp,'color','m','marker','o','linestyle','none'); end
	if any(iPn), hEp=line(TEMP(iPn,J(j)),PRES(iPn,J(j))); set(hEp,'color','g','marker','o','linestyle','none'); end
	axis ij; xlabel temp; ylabel pres; title(['CYCLE_NUMBER ',int2str(CYCLE_NUMBER(J(j)))],'interpreter','none');
	% Salinity plot
	aS=subplot(M,2,2*((j-ju)-1)+2); 
	hrS=plot(sal,pres,'color',grey); 
	hS=line(PSAL(:,jj),PRES(:,jj)); 
	heS=line(PSAL(:,J(j)),PRES(:,J(j))); 
	if any(iSo), hEp=line(PSAL(iSo,J(j)),PRES(iSo,J(j))); set(hEp,'color','m','marker','s','linestyle','none'); end
	if any(iPo), hEp=line(PSAL(iPo,J(j)),PRES(iPo,J(j))); set(hEp,'color','g','marker','s','linestyle','none'); end
	if any(iSn), hEp=line(PSAL(iSn,J(j)),PRES(iSn,J(j))); set(hEp,'color','m','marker','o','linestyle','none'); end
	if any(iPn), hEp=line(PSAL(iPn,J(j)),PRES(iPn,J(j))); set(hEp,'color','g','marker','o','linestyle','none'); end
	axis ij; xlabel sal; ylabel pres; 
	% Cosmetics:
	set([hT;hS],'color','b'); 
	set([heT;heS],'color','r','linestyle','-');
	% Monthly color and thickness gradient in refdata based on time difference in year:
	for i=2:11, dmon==i; set([hrT(ans),hrS(ans)],'color',brighten(grey,-.2+.1*i),'linewidth',1.5-.12*i); end
	for i=0:1,  dmon==i; set([hrT(ans),hrS(ans)],'color',brighten(grey+[0 .1 0],-.2+.1*i),'linewidth',1.5-.12*i); end
	% Inlay with positions:
	get(aT,'position'); ax=axes('position',[ans(1)+ans(3)-ans(3)/3-.01 ans(2)+.01 ans(3)/3 ans(4)/3]); hold on
	set(ax,'visible','off'); 
	plot(LO{j},LA{j},'color',grey);
	text(nanmean(LO{j}),nanmax(LA{j}),dtnote,'color',grey,'verticalalignment','bottom','horizontalalignment','center');
	plot(LONG(jj),LAT(jj),'marker','.','color','b'); 
	plot(LONG(J(j)),LAT(J(j)),'marker','s','color','r'); 
	% Print if page is full:
	if ismember(j,[5:5:N N]), 
	  print(gcf,'-depsc',[outfiles{I}(1:end-4),'_RTQCcheck',int2str(ju/5+1),'.eps']); 
	end 
  	if j==min(5,N) % If multipanel plots are made, generate part for report with just the first graph
	  fprintf(fid,'%s\n','Examples are shown in Figure~\ref{fig:RTQCcheck}. ');
	  fprintf(fid,'%s\n','\begin{figure}[H]');
	  fprintf(fid,'%s%u%s%u%s\n','\centerline{\includegraphics[width=1.2\textwidth,natwidth=',shape(3),',natheight=',shape(4),']{\floatsource\WMOnum_RTQCcheck1}}');
	  fprintf(fid,'%s\n','\caption{Float~\WMOnum. Profiles that have bad RTQC or DMQC flags (red) compared to float profiles before and after (blue) and nearby reference data (grey with darker shading and thicker lines for similar time of year, and green tint for same season).  Temperature in left and salinity in right panels. Flagged data are marked with magenta squares (RTQC) and circles (DMQC), while pressure flags have green marks. The inlays in left panels show the piece of float track in blue and the areas of origin of reference data as grey enclosures. These panels show the first flagged profiles for this float (or all; see Table~\protect\ref{tab:rtqc}). For geography, refer to Figure~\protect\ref{fig:float-info} Panel~1.}');
	  fprintf(fid,'%s\n','\label{fig:RTQCcheck}');
	  fprintf(fid,'%s\n','\end{figure}');
	end
      end % loop inspect or plot all profiles with flags
    else % No bad flags, then just a sentence for the report
      fprintf(fid,'%s\n','There were no RTQC flags to inspect nor any new flags made in DMQC for Float~\WMOnum.');
      system(['rm -f ',outfiles{I}(1:end-4),'_RTQCcheck*.eps']); % Clean up any previous plot
    end % Any bad flags
    fclose(fid);
  end % if new profiles
  
  % Check if there were more profiles than already checked:
  if checked{I}<n
    disp(['prepare_floats: New profiles downloaded. Update the parameter ''',checked_name,''' in INIT_DMQC from ',int2str(checked{I}),' to ',int2str(visd(2)),' for Float ',float_names{I},'!']);
    snippet('This float is still active and further monitoring is required.','monitoring');
  elseif direction=='A' & DATES(end)<dyear(now)-1/52*3 % three weeks allowance 
    ['There are no new profiles from this float since ',datestr(datenum(DATES(end),1,1),1),'.'];
    snippet(ans,'monitoring');
    disp(['prepare_floats: ',ans,' Has Float ',float_names{I},' been decommissioned after Cycle ',int2str(CYCLE_NUMBER(n)),' (profile ',int2str(n),')?']);
  end 
  clear LO LA
  
  % whos *qcn
  
  
  if direction=='D'
    % Re-insert the removed far-off surface values to the _ADJUSTED:
    PSAL_ADJUSTED(igi)=igiS; 
    % The application of bad flags comes below.
    % Some checks on this are: 
    % Data NaNed in the end and flag changes done:
    %[igi,PSAL_ADJUSTED(igi),str2num([PSALqco(igi),PSALqcn(igi)])]
    % Are all data NaN when flag is set to '4'?
    %[isnan(PSAL_ADJUSTED(igi)),str2num([PSALqco(igi),PSALqcn(igi)])]
  end

  
  
  
  %%%%%%% Time-series plots of the variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % PLOT the Hov-Möller diagrams of PTEMP and PSAL with raw data (i.e., disregarding flags)
  figure(1001); set(gcf,'OuterPosition',[1 385 1026 600]);
  ccc=1:max(CYCLE_NUMBER);								% All possible/expected cycle numbers
  [ppp,ttt,sss,ppq,ttq,ssq]=deal(nan(size(PRES,1),length(ccc)));			% Corresponding matrices (empty)
  ppp(:,CYCLE_NUMBER)=PRES; ttt(:,CYCLE_NUMBER)=PTMP; sss(:,CYCLE_NUMBER)=PSAL;		% Put data in right columns
  [ppq,ttq,ssq]=deal(repmat(' ',size(PRES,1),length(ccc)));				% Corresponding qc matrices (empty)
  ppq(:,CYCLE_NUMBER)=PRESqcn; ttq(:,CYCLE_NUMBER)=TEMPqcn; ssq(:,CYCLE_NUMBER)=PSALqcn;% Put data in right columns
  [~,X,Y]=profcolor(ccc,ppp,ttt);shading flat;						% Hov-Moller plot
  X=X(:,CYCLE_NUMBER); Y=Y(:,CYCLE_NUMBER);		% put coordinates back into compact format (for later marks on plot)
  axis ij; xlabel('CYCLE NUMBER'); ylabel('PRES');
  cmocean thermal; cbh=colorbar; ylabel(cbh,'PTMP');
  title(['Float ',char(float_names{I}),' Potential Temperature']);
  caxis(mima(ttt(ppp>10)));
  figure(1002); set(gcf,'OuterPosition',[1 385 1026 600]);
  profcolor(ccc,ppp,sss);shading flat; 
  axis ij; xlabel('CYCLE NUMBER'); ylabel('PRES')
  cmocean haline; cbh=colorbar; ylabel(cbh,'PSAL');
  title(['Float ',char(float_names{I}),' Salinity (PSS-78)']);
  caxis(mima(sss(ppp>10)));
  % Revisit these and print them in the end.
  
  % PLOT and PRINT the profile first-pressure series:
  figure(1003); set(gcf,'OuterPosition',[1 185 1026 900]);%[1 385 1026 900]);
  ap(1)=subplot(2,1,1); set(ap(1),'position',get(ap(1),'position')+[0 .15 0 -.15]);
  hP1=plot(CYCLE_NUMBER,PRES_ADJUSTED(1,:),'b.'); grid; axis ij % ylabel('Surface pressure (PRES_0)  [dbar]');
  PRESqco(1,:)=='2'; h=line(CYCLE_NUMBER(ans),PRES_ADJUSTED(1,ans));set(h,'color','c','marker','s','linestyle','none');
  PRESqco(1,:)=='3'; h=line(CYCLE_NUMBER(ans),PRES_ADJUSTED(1,ans));set(h,'color','m','marker','s','linestyle','none');
  PRESqco(1,:)=='4'; h=line(CYCLE_NUMBER(ans),PRES_ADJUSTED(1,ans));set(h,'color','r','marker','s','linestyle','none');
  ap(2)=subplot(2,1,2); set(ap(2),'position',get(ap(2),'position')+[0 0 0 .15+.1]);
  CYC=repmat(CYCLE_NUMBER,size(PRESqco,1),1);
  [1:10:m m]; PRESqco10=PRESqco(ans,:); PRES10=PRES_ADJUSTED(ans,:); CYC10=CYC(ans,:);
  hP2=plot(CYC10',PRES10','b.'); grid; axis ij
  set(ap(1),'XAxisLocation','top');
  xlabel(ap(2),'CYCLE NUMBER'); ylabel(ap,'PRES  [dbar]');
  xlim(ap,[.5 max(CYCLE_NUMBER)+.5]);
  get(ap(1),'ylim'); set(ap(1),'ylim',[min(-1,ans(1)-1) max(5,ans(2)+1)]);
  get(ap(2),'ylim'); set(ap(2),'ylim',[min(-1,ans(1)-1) profdepth+500]);
  PRESqco10=='2'; hPo{2}=line(CYC10(ans)',PRES10(ans)'); set(hPo{2},'color','c','marker','s','linestyle','none','tag','qcflag','userdata','RTQC flag ''2''');
  PRESqco10=='3'; hPo{3}=line(CYC10(ans)',PRES10(ans)'); set(hPo{3},'color','m','marker','s','linestyle','none','tag','qcflag','userdata','RTQC flag ''3''');
  PRESqco=='4'; hPo{4}=line(CYC(ans)',PRES_ADJUSTED(ans)'); set(hPo{4},'color','r','marker','s','linestyle','none','tag','qcflag','userdata','RTQC flag ''4''');
  PRESqco=='8'; hPo{8}=line(CYC(ans)',PRES_ADJUSTED(ans)'); set(hPo{8},'color','b','marker','s','linestyle','none','tag','qcflag','userdata','RTQC flag ''8''');
  %PRESqco=='9'; hPo{9}=line(CYC(ans)',PRES_ADJUSTED(ans)'); set(hPo{9},'color','k','marker','s','linestyle','none','tag','qcflag','userdata','RTQC flag ''9''');
  % Flag '9' is for missing data, so will not show anyway.
  title(ap(1),['Top of profile pressure measurements by float ',char(float_names{I})]);
  title(ap(2),['Every 10th pressure measurement by float ',char(float_names{I})]);
  % Legend and print of lower panel done below, after DMQC.
  
  
        % ----- APPLY FLAGS TO THE DATA: --------------------------------
  
	% Apply flags from the RTQC to the _ADJUSTED: 
	% For reversal qcn has been set to 1. Do not reverse the RTQC flags
	% yet, because they should be shown properly in the graphics and
	% stats, but make sure reversed flags are not applied to the data here
	% (added '&<PARAMETER>qcn~=1' to this code).
	% NaN out columns based on bad coordinates:
	%j=find(POSqco=='4'&POSqcn~='1' | JULDqco=='4'&JULDqcn~='1' & LAT<-90 & 90<LAT);  
	%j=find(POSqco=='4'&POSqcn~='1' | JULDqco=='4'&JULDqcn~='1' | LAT<-90 | 90<LAT);  
	j=find(POSqco=='4' & POSqcn~='1' | JULDqco=='4' & JULDqcn~='1');  
	if any(j), 
	  disp(['prepare_floats: Direction ''',direction,''' profiles from cycles ',zipnumstr(CYCLE_NUMBER(j)),...
		' have been NaN''ed due to RTQC (unreversed)  bad date or position!'])
	  LONG(j)=NaN; LAT(j)=NaN; DATES(j)=NaN; %PRES_ADJUSTED(:,j)=NaN; PSAL_ADJUSTED(:,j)=NaN; TEMP_ADJUSTED(:,j)=NaN;
	end
	% NaN out single values (according to manual Section 3.3):
	PRESqco=='4' & PRESqcn~='1'; PRES_ADJUSTED(ans)=NaN; 
	PSAL_ADJUSTED(ans | PSALqco=='4' & PSALqcn~='1')=NaN; 
	TEMP_ADJUSTED(ans | TEMPqco=='4' & TEMPqcn~='1')=NaN;
	% Also remove the ASD profiles prior to the comparison and final
        % status plots. (This does not affect the flag-objects, hense
        % not the statistics or plotmarks, and data will be re-inserted
        % before OWC.):
  	if ASD<=n, PSAL_ADJUSTED(:,ASD:end)=NaN; end

	% Apply flags from the DMQC to the _ADJUSTED: 
	% Now apply new bad-flags to the data going into OWC, i.e. _ADJUSTED.
	% This is where we need to heed the fact that qcn now contains
        % _only_the_brand_new_ flags.
	% Reversed data is given qcn=1 and considered above when removing
	% data based on RTQC. Now we only need to NaN out the new '4's.
	% NaN out columns:
	%j=find(POSqcn=='4' | JULDqcn=='4' & LAT<-90 & 90<LAT); 
	%j=find(POSqcn=='4' | JULDqcn=='4' | LAT<-90 | 90<LAT); 
	j=find(POSqcn=='4' | JULDqcn=='4'); 
	if any(j), 
	  disp(['prepare_floats: Direction ''',direction,''' profiles from cycles ',zipnumstr(CYCLE_NUMBER(j)),...
		' have been NaN''ed due to bad date or position found in DMQC!'])
	  LONG(j)=NaN; LAT(j)=NaN; DATES(j)=NaN; %PRES_ADJUSTED(:,j)=NaN; PSAL_ADJUSTED(:,j)=NaN; TEMP_ADJUSTED(:,j)=NaN;
	end
	% NaN out single values (according to manual, Section 3.3):
	PRESqcn=='4'; PRES_ADJUSTED(ans)=NaN; PSAL_ADJUSTED(ans | PSALqcn=='4')=NaN;  TEMP_ADJUSTED(ans | TEMPqcn=='4')=NaN;
	% At this point the variables here are to be considered _ADJUSTED variables.
	% PRES, TEMP, and PSAL are never changed in this script.

	% whos LONG LAT DATES *_ADJUSTED

	
	
  %%%%%	Phase 5: Overall comparison with the reference data %%%%%%
  
  [yrnow,~]=datevec(now);
  long=LONG; long>180; long(ans)=long(ans)-360;
  m_proj('albers','lon',mima(long)+[-1 1],'lat',mima(LAT)+[-.5 .5]);
  [x,y]=m_ll2xy(long,LAT);
  % Find and put cluster numbers (i.e., groups) in sequence:
  if isscalar(x),	CLU1=1; Nclu{I}=1;					% Only one cluster possible
    msg='prepare_floats: Only one cluster possible';
  elseif ischar(Nclu{I})							% Use the calseries to cluster
    if direction=='A' % & exist('calseries','var')
      msg=['prepare_floats: Using calseries to define clusters. '];
      CLU1=calseries;
      Nclu{I}=length(unique(CLU1));
    else
      msg=['prepare_floats: The use of calseries to define clusters are only valid for ascending profiles. '];
      let2int(Nclu{I});
      if any(ans)						% Number given by letter
	msg=[msg,'Using ',int2str(ans(1)),' clusters as given by the character ''',Nclu{I},''' in Nclu.'];
	Nclu{I}=ans(1); CLU1=kmeans([x',y'],Nclu{I})';		
      else
	CLU1=1; Nclu{I}=1;					% Default for descending profiles
	msg=[msg,'Defaulting to one cluster.'];
      end
    end
    return
  elseif iscell(Nclu{I})							% Cluster manually
    msg=['prepare_floats: Using manual clustering. '];
    % cell2mat(Nclu{I}) is a vector of split-points
    CLU1=nan(1,n);
    [0 cell2mat(Nclu{I}) n]; for i=1:length(ans)-1, CLU1(ans(i)+1:ans(i+1))=i; end
    Nclu{I}=length(cell2mat(Nclu{I}))+1;
  elseif imag(Nclu{I}),	
    msg=['prepare_floats: Automatic clustering wrt. time. '];
    CLU1=kmeans([CYCLE_NUMBER',x',y'],abs(Nclu{I}))';	% Cluster wrt time.
  else
    msg=['prepare_floats: Clustering is automatic. '];
    CLU1=kmeans([x',y'],abs(Nclu{I}))';			% Cluster not wrt time.
  end
  disp(msg);
  CL=unique(CLU1,'stable'); CL=CL(~isnan(CL)); 
  CLU=nan(size(CLU1)); 
  for j=1:length(CL), CLU(CLU1==CL(j))=j; end 
  % Prepare the two figure windows
  figure(1004); figure(1005); shape=[1 1 800 250*abs(Nclu{I})];
  set([1004],'position',shape,'paperunits','points','paperposition',shape,'PaperOrientation','portrait','PaperPositionMode','manual');
  set([1005],'position',shape+[800 0 0 0],'paperunits','points','paperposition',shape,'PaperOrientation','portrait','PaperPositionMode','manual');
  clear a h ax ja
  for j=1:abs(Nclu{I})
    jj=find(CLU==j); 
    %if isempty(jj), continue; end
    tit=['Cluster ',int2str(j),' (Cycles ',zipnumstr(CYCLE_NUMBER(jj)),' ; months ',zipnumstr(MONTHS(jj)),')'];
    [LO{j},LA{j}]=halo(LONG(jj),LAT(jj),.3,.2); 
    [pres,temp,sal,long,lat,dates,~,dt,dmon,RMONTHS]=inpolygon_referencedata(LO{j},LA{j},tardir,0,time(jj));
    save(['refdata_cluster',int2str(j)],'pres','temp','sal','long','lat','dates');
    figure(1004); % Profile plots
    plim=[MAP_P_EXCLUDE max(PRES_ADJUSTED(:,jj),[],'all')];
    if diff(plim)<50 % In case of shallow profiles and/or ow_config.txt not updated properly
      warning(['MAP_P_EXCLUDE ',sprintf('%d',round(plim(1))),...
	       ' too large for profiles of max pressure ',sprintf('%d',round(plim(2))),...
	       ' dbar. Please edit ow_config.txt for float ',float_names{I},'.']);
      plim(1)=0; 
    end 
    a(j,1)=subplot(abs(Nclu{I}),2,2*(j-1)+1); 
    hrT=plot(temp,pres,'color',grey); 
    for i=0:5			% Monthly color gradient in refdata based on time difference in year
      dcm=brighten(grey,-.2+.1*i); dt>=i*30; set(hrT(ans),'color',dcm); % Darker near, brighter far 
    end
    h=line(TEMP_ADJUSTED(:,jj),PRES_ADJUSTED(:,jj)); set(h,'color','b');
    %get(a(j,1),'ylim'); set(a(j,1),'ylim',[MAP_P_EXCLUDE ans(2)]);
    set(a(j,1),'ylim',plim);
    axis ij; xlabel temp; ylabel pres; title(tit);
    a(j,2)=subplot(abs(Nclu{I}),2,2*(j-1)+2); 
    hrS=plot(sal,pres,'color',grey);
    %get(a(j,2),'ylim'); set(a(j,2),'ylim',[MAP_P_EXCLUDE ans(2)]);
    for i=0:5			% Monthly color gradient in refdata based on time difference in year
      dcm=brighten(grey,-.2+.1*i); dt>=i*30; set(hrS(ans),'color',dcm); % Darker near, brighter far 
    end
    h=line(PSAL_ADJUSTED(:,jj),PRES_ADJUSTED(:,jj)); set(h,'color','b');
    set(a(j,2),'ylim',plim);
    axis ij; xlabel sal; ylabel pres; title(tit);
    %
    % Informative text about refdata in cluster:
    if ~isempty(pres)
      dtnote={'Reference data are from ';['months ',zipnumstr(RMONTHS),' ']};
      % Mean refdata salinity at cold, deep ptemp layers:
      ptmp=sw_ptmp(sal,temp,pres,0);
      indt=find(ptmp<-0.7 & pres>1000);
      if ~isempty(indt)
	Msal=round(nanmean(sal(indt)),3,'decimals'); 
	Mpres=round(min(pres(indt))); 
	dtnote{3}=['Mean refsal = ',num2str(Msal,'%6.3f'),' @ pres > ',int2str(Mpres),' db '];
      end
      text(max(xlim),max(ylim),dtnote,'color','k','verticalalignment','bottom','horizontalalignment','right');
    end
    %
    % Trend plots (all year ref data used here):
    figure(1005);  
    [pres,temp,sal,long,lat,dates]=inpolygon_referencedata(LO{j},LA{j},tardir);
    pres<MAP_P_EXCLUDE; pres(ans)=NaN; temp(ans)=NaN; sal(ans)=NaN;  % as these depths are ignored in the OWC-analysis
    tit=['Cluster ',int2str(j),' (Cycles ',zipnumstr(CYCLE_NUMBER(jj)),' ; all year)'];
    ja(j,1)=subplot(abs(Nclu{I}),2,2*(j-1)+1); 
    [years,~]=datevec(rdtime(dates)); years=years(:)';
    g_y=yrnow-30:yrnow+1; % center around newyear/winter 
    % Flexible depth grid:
    %%g_d=200:400:2000; 
    g_d=MAP_P_EXCLUDE+200:400:2000; 
    find(min(PRES_ADJUSTED,[],'all') < g_d & g_d < max(PRES_ADJUSTED,[],'all'));
    if length(ans)==1, ans(1)+[0 1]; elseif isempty(ans), ans=1:2; end
    g_d=g_d(ans);
    %if max(pres,[],'all')<800,	g_d=[200 600]; 
    %else			g_d=1000:400:2000; end
    bin=bin2d(repmat(years(:)',size(sal,1),1),pres,sal,g_y,g_d,false);
    jbin=bin2d(repmat(DATES(jj),size(PSAL_ADJUSTED(:,jj),1),1),PRES_ADJUSTED(:,jj),PSAL_ADJUSTED(:,jj),g_y,g_d,false);
    %errorbar(repmat(bin.x,size(bin.mean,1),1)',bin.mean',bin.var'./(bin.n'-1)); % error of the mean 
    errorbar(repmat(bin.x,size(bin.mean,1),1)',bin.mean',bin.var'); % variance 
    h=line(jbin.x,jbin.mean); set(h,'linestyle','none','marker','*')
    title(tit); ylabel sal; grid; xlim([g_y(1)-1,g_y(end)+1]);
    if j==1,al=legend([strcat({'Reference data '},int2str(jbin.yg(1:end-1)'),'-',int2str(jbin.yg(2:end)'),{' m'});...
		       strcat({'Float data '},int2str(jbin.yg(1:end-1)'),'-',int2str(jbin.yg(2:end)'),{' m'})]); 
    end 
    if all(isnan(jbin.mean),'all')|all(isnan(bin.mean),'all') 
      text(mean(xlim),mean(ylim),'No data to compare!','HorizontalAlignment','center','BackgroundColor',[.9 .9 .9]);
      if j==1, set(al,'visible','off'); end
    else
      %mima(bin.mean'+bin.var'./(bin.n'-1),bin.mean'-bin.var'./(bin.n'-1),jbin.mean); % error of the mean
      mima(bin.mean'+bin.var',bin.mean'-bin.var',jbin.mean); % variance
      if diff(ans)<0.01, ans+[-.005 .005]; end
      ylim(ans);
    end
  end
  % Put cluster-map in plots:
  figure(1004);
  get(a(1),'position'); ax=axes('position',[ans(1)+ans(3)-ans(3)/3-.01 ans(2)+.01 ans(3)/3 ans(4)/3]); hold on
  set(ax,'visible','off');
  for j=1:abs(Nclu{I})
    plot(LO{j},LA{j},'color',grey); text(LO{j}(1),LA{j}(1),int2str(j),'color',grey);
  end
  plot(LONG,LAT,'marker','.','color','b'); 
  ax2=copyobj(ax,1005);							% Cluster map also on trend plot
  get(al,'position'); set(al,'position',[0.8-ans(3) ans(2:4)]);		% Move legend out of the way
  get(ax2,'position'); set(ax2,'position',[0.6-ans(3) ans(2:4)]);	% Move cluster map out of the way
  round(mima(get(ja,'ylim'))*100)/100; %set(ja,'ylim',ans);
  set(ja,'ytick',ans(1):.01:ans(2),'yminortick','on','yminorgrid','on');
  % Print figures to files:
  print(1004,'-depsc',[outfiles{I}(1:end-4),'_refcomp.eps'])
  print(1005,'-depsc',[outfiles{I}(1:end-4),'_trendcheck.eps'])

  % Generate figure parts for report:
  fid=fopen('refcomp.tex','w');
  fprintf(fid,'%s%u%s%u%s\n','\centerline{\includegraphics[width=1.2\textwidth,natwidth=',shape(3),',natheight=',shape(4),']{\floatsource\WMOnum_refcomp}}');
  fprintf(fid,'%s\n',['\caption{Float \WMOnum\ compared to nearby reference data. Float profiles (blue lines) are divided into clusters based on positions, and compared to nearby profiles from the reference data set (grey lines with darker shading for similar time of year). Temperature in left and salinity in right panels, and one row per cluster. The upper ',int2str(MAP_P_EXCLUDE),'~m are omitted as these depths are ignored in the OWC-analysis. The inlay in first panel shows areas of origin of reference data used as grey enclosures around the clusters of positions on the blue float track.  For geography, refer to Figure~\protect\ref{fig:float-info} Panel~1.  Note that for floats with ASD, uncorrectable salinity profiles are not shown here. }']);

  fprintf(fid,'%s\n','\label{fig:refcomp}');
  %fprintf(fid,'%s\n','\end{figure}');
  fclose(fid);
  fid=fopen('trendcheck.tex','w');
  %fprintf(fid,'%s\n','\begin{figure}[H]');
  %fprintf(fid,'%s\n','\centering');
  fprintf(fid,'%s%u%s%u%s\n','\centerline{\includegraphics[width=1\textwidth,natwidth=',shape(3),',natheight=',shape(4),']{\floatsource\WMOnum_trendcheck}}');
  fprintf(fid,'%s\n','\caption{Temporal evolution of data from Float~\WMOnum\ and reference data (by cluster as in Figure~\protect\ref{fig:refcomp}). Time series of annual means of reference data in depth bins (see legend) are plotted as coloured lines with error bars representing the bin variance, and annual bin means of float data by the same method plotted with asterisks. Annual bins are centered around new year (i.e., winter). The sketch under the legend shows areas of origin of reference data used as grey enclosures around the clusters of positions on the blue float track.  For geography, refer to Figure~\protect\ref{fig:float-info} Panel~1.  }');
  fprintf(fid,'%s\n','\label{fig:trendcheck}');
  %fprintf(fid,'%s\n','\end{figure}');
  fclose(fid);
  clear LO LA
  
 

  

  %%%%%%%%%%%%%%%% CORRECT DEEP ARGO/ARVOR FLOATS PRESSURE DEPENDENT CONDUCTIVITY BIAS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % "CPcor_new, should be used to re-compute salinity as the first step in 'D' mode."
  fid=fopen('prescondbias.tex','w');
  if strcmp(platyp,'ARVOR-D') & logical(1) & direction=='A'
    disp(['prepare_floats: Adjusting pressure dependent conductivity bias in ',ctdmodel,' on ',platyp,' float.']);
    % First create dummy raw parameters affected by flags, for use here: 
    PSALf=PSAL_ADJUSTED; TEMPf=TEMP_ADJUSTED; PRESf=PRES_ADJUSTED; 
	
    % 1. Use a new CPcor value, CPcor_new, to re-compute salinity
    % The steps to re-compute salinity with CPcor_new are as follows.
    % (a). Fill PRES_ADJUSTED and TEMP_ADJUSTED.
    % Follow the same 'D' mode procedures as for the 2000-dbar floats,
    % as in Sections 3.3 and 3.4.    (det er før OWC)
    % √ Filled above already.
    % (b). Compute original conductivity, Co.  This is done by using
    % PRES, TEMP, PSAL. For example, if using the Gibbs-SeaWater (GSW)
    % Oceanographic Matlab Toolbox, then Co = gsw_C_from_SP (PSAL, TEMP,
    % PRES).
    ii=find(~isnan(PRESf)&~isnan(PSALf)&~isnan(TEMPf)&~isnan(PRES_ADJUSTED)&~isnan(TEMP_ADJUSTED));	% All parameters non-NaN because 
    %Co=gsw_C_from_SP(PSAL(ii),TEMP(ii),PRES(ii));	% gsw_ can't take any NaNs! :-(
    Coo=gsw_C_from_SP(PSALf(ii),TEMPf(ii),PRESf(ii));	% gsw_ can't take any NaNs! :-(
    Co=nan(m,n); Co(ii)=Coo; clear Coo 		% Matrix again
    % (c). Compute new conductivity, Cnew.
    % Cnew = Co*(1 + d*TEMP + CPcor_SBE*PRES) / (1 + d*TEMP_ADJUSTED + CPcor_new*PRES_ADJUSTED); 
    % where d = 3.25e-06,
    d = 3.25e-6;
    % CPcor_SBE = -9.57e-08 dbar-1,
    % CPcor_new = -12.5e-08 dbar-1 for SBE-61 data,
    % CPcor_new = -13.5e-08 dbar-1 for Deep SBE-41CP data, or
    % CPcor_new = as estimated by DMQC operator.
    CPcor_SBE = -9.57e-8;
      switch ctdmodel
       case 'SBE41CP', CPcor_rec = -13.5e-8; % recommended
       case 'SBE61',  CPcor_rec = -12.5e-8;
      end
    try
      load(['CPcor',filesep,'operatorCPcor']); % Load estimated CPcor_new
      % about=['For this float the operator has determined CPcor\_new and cell-gain based on a near-deployment ship CTD profile. ',...
      % 	     'The cell-gain is a factor on the conductivity that also removes offset relative to the CTD-profile.'];
      about='For this float the operator has determined CPcor\_new based on a near-deployment ship CTD profile. ';
      endtext=' compared to a near-deployment ship CTD profile. ';
      disp(['prepare_floats: ',about]);
      capt=['Float \WMOnum. Salinity profiles from the ',num2ordinal(CYCLE_NUMBER(np)),' cycle',...
	    ' using the manufacturer calibration ($CPcor_{SBE} =$ ',num2str(CPcor_SBE*1.0e+8,'%6.2f'),'~$\cdot 10^{-8}~dbar^{-1}$; i.e., raw data; green),',...
	    ' using the recommended compressibility term ($CPcor_{new} =$ ',num2str(CPcor_rec*1.0e+8,'%6.2f'),'~$\pm1.5~\cdot~10^{-8}~dbar^{-1}$ ; cyan$\pm$dashed),',...
	    ' and using the operator''s compressibility term ($CPcor_{new} =$ ',num2str(CPcor_new*1.0e+8,'%6.2f'),'~$\cdot~10^{-8}~dbar^{-1}$ ; magenta).',...
	    ' The ship CTD profile used is from Cruise ',ctd_cruise,', Station ',ctd_station,' (black).',...
	    ' Dashed horisontal line shows the minimum pressure the fit is done for.'];

      % Calibration info (Note: "please record the CPcor correction first, then other delayed-mode salinity steps next."):
      scientific_calib.equation.PSAL	= repmat(['New conductivity = original conductivity * (1 + delta*TEMP + CPcor_SBE*PRES) / (1 + delta*TEMP_ADJUSTED + CPcor_new*PRES_ADJUSTED). '],n,1); % To be placed first!
      scientific_calib.coefficient.PSAL	= repmat(['CPcor_new = ',sprintf('%0.3g',CPcor_new),'; CPcor_SBE = ',sprintf('%0.3g',CPcor_SBE),'; delta = ',sprintf('%0.3g',d),'. '],n,1); % To be placed first!
      scientific_calib.comment.PSAL	= repmat(['New conductivity computed by using a different CPcor value from that provided by Sea-Bird. '],n,1); % To be placed first!
    
    catch
      disp('prepare_floats: No (valid) operator CPcor_new found.');  
      % Should have found 'ctd_*','CPcor_new','leg','np','minPRESS');
      find(max(PRES_ADJUSTED)>2500); np=ans(1); % Find first profile deep enough
      leg.float=[datestr(time(np)),', ',num2str(LONG(np),'%10.3f'),'\circE, ',num2str(LAT(np),'%10.3f'),'\circN, Argo ',float_names{I},', cycle ',int2str(np)]
      leg.ctd='no CTD';
      minPRESS=nan;
      % ctd_* are not used below in the case of no operator CPcor_new.
      CPcor_new=999; % so the below IF test fails
    end
    %     -20e-08 dbar -1 <= CPcor_new <= - 5e-08 dbar -1
    % If CPcor_new falls outside those ranges, the DM operator should consider that the salinity
    % data are not correctable and that further investigations are required to demonstrate that
    % such CPcor_new values are valid.
    if CPcor_new <= -20e-8 | -5e-8 <= CPcor_new
      CPcor_new=CPcor_rec % use the recommended
      about=['For this float the operator has not determined, or not been able to determine, any valid CPcor\_new. ',...
	    'Instead the recommended CPcor\_new for ',ctdmodel,' is used.'];
      endtext='. ';
      disp(['prepare_floats: ',about]);
      capt=['Float \WMOnum. Salinity profiles from the ',num2ordinal(CYCLE_NUMBER(np)),' cycle',...
	    ' using the manufacturer calibration ($CPcor_{SBE} =$ ',num2str(CPcor_SBE*1.0e+8,'%6.2f'),'~$\cdot 10^{-8}~dbar^{-1}$; i.e., raw data; green),',...
	    ' and using the recommended compressibility term ($CPcor_{new} =$ ',num2str(CPcor_rec*1.0e+8,'%6.2f'),'~$\pm1.5~\cdot~10^{-8}~dbar^{-1}$ ; cyan$\pm$dashed).',...
	    ' For this float the recommended term is used (magenta).'];

      % Calibration info:
      scientific_calib.equation.PSAL	= repmat(['New conductivity = original conductivity * (1 + delta*TEMP + CPcor_SBE*PRES) / (1 + delta*TEMP_ADJUSTED + CPcor_new*PRES_ADJUSTED). '],n,1); % To be placed first!
      scientific_calib.coefficient.PSAL	= repmat(['CPcor_new = ',sprintf('%0.3g',CPcor_new),'; CPcor_SBE = ',sprintf('%0.3g',CPcor_SBE),'; delta = ',sprintf('%0.3g',d),'. '],n,1); % To be placed first!
      scientific_calib.comment.PSAL	= repmat(['New conductivity computed by using a different CPcor value from that provided by Sea-Bird. '],n,1); % To be placed first!
    
    end
    %
    % (c). Compute new conductivity, Cnew.
    Cnew = Co .* (1 + d*TEMP + CPcor_SBE*PRES) ./ (1 + d*TEMP_ADJUSTED + CPcor_new*PRES_ADJUSTED); 
    % Show recommended also (e):
    Cnew_recu = Co .* (1 + d*TEMP + CPcor_SBE*PRES) ./ (1 + d*TEMP_ADJUSTED + (CPcor_rec+1.5e-8)*PRES_ADJUSTED); 
    Cnew_rec  = Co .* (1 + d*TEMP + CPcor_SBE*PRES) ./ (1 + d*TEMP_ADJUSTED + CPcor_rec*PRES_ADJUSTED); 
    Cnew_recl = Co .* (1 + d*TEMP + CPcor_SBE*PRES) ./ (1 + d*TEMP_ADJUSTED + (CPcor_rec-1.5e-8)*PRES_ADJUSTED); 
    %Cnew = M_new * Co .* (1 + d*TEMP + CPcor_SBE*PRES) ./ (1 + d*TEMP_ADJUSTED + CPcor_new*PRES_ADJUSTED);  % (This with cell-gain M_new. M_new=1 when no operator CPcor_new is made.)
    % (d). Compute new salinity, PSAL_ADJUSTED_Cnew.  This is done by
    % using PRES_ADJUSTED, TEMP_ADJUSTED, and Cnew. For example, if
    % using the Gibbs-SeaWater (GSW) Oceanographic Matlab Toolbox, then
    % PSAL_ADJUSTED_Cnew = gsw_SP_from_C (Cnew, TEMP_ADJUSTED, PRES_ADJUSTED).
    PSAL_ADJUSTED_Cnew=nan(size(PSAL_ADJUSTED));
    PSAL_ADJUSTED_Cnew(ii)=gsw_SP_from_C(Cnew(ii),TEMP_ADJUSTED(ii),PRES_ADJUSTED(ii));
    % Show recommended also (e):
    [PSAL_ADJUSTED_Cnew_recu,PSAL_ADJUSTED_Cnew_rec,PSAL_ADJUSTED_Cnew_recl]=deal(nan(size(PSAL_ADJUSTED)));
    PSAL_ADJUSTED_Cnew_recu(ii)=gsw_SP_from_C(Cnew_recu(ii),TEMP_ADJUSTED(ii),PRES_ADJUSTED(ii));
    PSAL_ADJUSTED_Cnew_rec(ii) =gsw_SP_from_C(Cnew_rec(ii), TEMP_ADJUSTED(ii),PRES_ADJUSTED(ii));
    PSAL_ADJUSTED_Cnew_recl(ii)=gsw_SP_from_C(Cnew_recl(ii),TEMP_ADJUSTED(ii),PRES_ADJUSTED(ii));
    %
    figure(2001);clf;set(gcf,'OuterPosition',[1027 385 650 700]);	% Look at correction of first profile.  
    dline=minPRESS*[1 1];
    leg.par.old=[' (CPcorSBE = ',num2str(CPcor_SBE*1.0e+8,'%6.2f'),' \cdot 10^{-8} dbar^{-1})'];
    leg.par.new=[' (CPcor\_new = ' num2str(CPcor_new*1.0e+8,'%6.2f') ' \cdot 10^{-8} dbar^{-1})'];
    leg.par.rec=[' (CPcor\_new = ' num2str(CPcor_rec*1.0e+8,'%6.2f') '\pm1.5 \cdot 10^{-8} dbar^{-1})'];
    hr=plot(PRES_ADJUSTED(:,np),PSAL_ADJUSTED(:,np),'g-',...
	    PRES_ADJUSTED(:,np),PSAL_ADJUSTED_Cnew_recu(:,np),'c--',...
	    PRES_ADJUSTED(:,np),PSAL_ADJUSTED_Cnew_rec(:,np),'c-',...
	    PRES_ADJUSTED(:,np),PSAL_ADJUSTED_Cnew_recl(:,np),'c--',...
	    PRES_ADJUSTED(:,np),PSAL_ADJUSTED_Cnew(:,np),'m-'); 
    clear *_rec*
    if strcmp(endtext,'. ') % no CTD-data
      legend([hr(1);hr(3);hr(5)],['PSAL ',leg.float,leg.par.old],...
	     ['PSAL\_ADJUSTED\_Cnew recommended for sensor model ',ctdmodel,leg.par.rec],...
	     ['PSAL\_ADJUSTED\_Cnew',leg.par.new],...
	     'location','north');
    else
      hc=line(ctd_z,ctd_salt); set(hc,'linestyle','-','color','k'); 
      hlim=line(dline,[33 36]); set(hlim,'linestyle','--','color','k'); % Line showing depth of CPcor fit.
      legend([hr(1);hc;hr(3);hr(5)],['PSAL ',leg.float,leg.par.old],...
	     ['CTD   ',leg.ctd],...
	     ['PSAL\_ADJUSTED\_Cnew recommended for sensor model ',ctdmodel,leg.par.rec],...
	     ['PSAL\_ADJUSTED\_Cnew',leg.par.new],...
	     'location','north');
    end
    ylim([34.895 34.925]); xlim([0 3500]);
    grid; xlabel PRES; ylabel PSAL;  view([90 90]);
    %title(['Pressure dependent conductivity-bias correction (Cycle ',int2str(CYCLE_NUMBER(np)),')']);
    title(['Pressure dependent conductivity-bias correction']);
    print(gcf,'-depsc',[outfiles{I}(1:end-4),'_prescondbias.eps']);
    
    % Create text-section for the report:
    %fprintf(fid,'%s\n','\newpage');
    fprintf(fid,'%s\n','\subsection{Correcting deep ARVOR float pressure dependent conductivity bias}');
    fprintf(fid,'%s%s%s%s%s%s%s\n','For ARVOR-D floats, conductivity has to be corrected due to a too small conductivity cell compressibility term (\emph{CPcor}) in the manufacturer calibration of the SBE CTDs. ',about,' Figure~\ref{fig:prescondbias} shows the uncorrected and corrected salinities from the ',num2ordinal(CYCLE_NUMBER(np)),' cycle',endtext,' If the first cycle is not used it is probably because it is too shallow (see Figure~\ref{fig:hovmoller}b).');
    %fprintf(fid,'%s\n','This re-computation salinity should be done as the first step in D-mode. However, to ensure the correct cycle is used for this calibration, it is done after the delayed-mode procedures for coordinates (Section~\ref{sec:DM-coordinates}). ');
    fprintf(fid,'%s\n','\begin{figure}[H]');
    fprintf(fid,'%s\n','\centering');
    fprintf(fid,'%s\n','\includegraphics[width=0.7\textwidth,natwidth=560,natheight=550]{\floatsource\WMOnum_prescondbias}');
    fprintf(fid,'%s%s%s\n','\caption{',capt,'}');
    fprintf(fid,'%s\n','\label{fig:prescondbias}');
    fprintf(fid,'%s\n','\end{figure}');
    
    % 2. Assess sensor drift and offset in 'D' mode by using PSAL_ADJUSTED_Cnew
    PSAL=PSAL_ADJUSTED;			% Save PSAL for the making of operator Cnew
					% and if no adjustment is being made.
    PSAL_ADJUSTED=PSAL_ADJUSTED_Cnew;	% For OWC and PSAL_ADJUSTED if no
					% further adjustment is done.
    
  else % not deep
    PSAL=PSAL_ADJUSTED;			% Save PSAL for the making of operator Cnew
					% and if no adjustment is being made.
    PSAL_ADJUSTED_Cnew=[];		% Empty for saving and write_D
  end % if deep arvor
  fclose(fid);

  % whos PSAL PSAL_ADJUSTED PSAL_ADJUSTED_Cnew
  
  
  % ----- Add SP adjustments to pressure figure: -------------
  if exist('SP') % PRES has been adjusted for surface pressure
    axes(ap(1));
    hSP(1)=line(SP.CN,SP.TPV); set(hSP(1),'marker','*','color','r')
    hSP(2)=line(SP.CN,SP.SP); set(hSP(2),'marker','o','color','r')
    set(hSP(1),'linestyle','none'); 
    legend([hP1 hSP],'Top of profile pressure (corrected)','Surface pressure offset','Despiked surface pressure offset','location','best');
    % Add PSAL scientific calib after potential CP_corr corrections:
    'Salinity re-calculated by using PRES_ADJUSTED and recorded in PSAL_ADJUSTED. ';
		scientific_calib.comment.PSAL(SP.IA,end+[1:size(ans,2)]) = repmat(ans,length(SP.IA),1);
  end
  
  
  
  %%%%%%% FINAL PHASE: Final plots, snippets, save for OWC, and alter R-files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  DENS = sw_pden(PSAL_ADJUSTED,TEMP_ADJUSTED,PRES_ADJUSTED,0);		% Update density for the final plot:
  
  plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'.eps']);	% Finally plot the clean plot.
  
  PTMP = sw_ptmp(PSAL_ADJUSTED,TEMP_ADJUSTED,PRES_ADJUSTED,0);		% Create this new variable

  %%%%%% Add flags to the Hov-Möller diagrams and the pressure plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % color:'w' variable itself
  % 	'grey' pressure
  % marker:point RTQC='2' or '3' 
  % 	square = RTQC='4'
  % 	diamond = reversed to '1'
  % 	circle = new '4'
  figure(1001); 
  if any(~isnan(PTMP),'all'), caxis(mima(PTMP)); end  % Change colourscale to the good data range
  find(TEMPqco=='2'); ho=line(X(ans),Y(ans)); set(ho,'marker','.','markersize',2,'color','w','linestyle','none');
  find(TEMPqco=='3'); ho=line(X(ans),Y(ans)); set(ho,'marker','.','markersize',2,'color','w','linestyle','none');
  find(TEMPqco=='4'); ho=line(X(ans),Y(ans)); set(ho,'marker','s','markersize',8,'color','w','linestyle','none');
  find(TEMPqcn=='1');   hr=line(X(ans),Y(ans)); set(hr,'marker','d','markersize',8,'color','w','linestyle','none');
  find(TEMPqcn=='4');   hn=line(X(ans),Y(ans)); set(hn,'marker','o','markersize',8,'color','w','linestyle','none');
  % find(PRESqco=='4'); ho=line(X(ans),Y(ans)); set(ho,'marker','s','markersize',8,'color',[.5 .5 .5],'linestyle','none');
  % find(PRESqcn=='1');   hr=line(X(ans),Y(ans)); set(hr,'marker','d','markersize',8,'color',[.5 .5 .5],'linestyle','none');
  % find(PRESqcn=='4');   hn=line(X(ans),Y(ans)); set(hn,'marker','o','markersize',8,'color',[.5 .5 .5],'linestyle','none');
  print(gcf,'-dpng',[outfiles{I}(1:end-4),'_HMtemp.png']);
  figure(1002); 
  if any(~isnan(PSAL_ADJUSTED),'all'), caxis(mima(PSAL_ADJUSTED)); end % Change colourscale to the good data range
  find(PSALqco=='2');  ho=line(X(ans),Y(ans)); set(ho,'marker','.','markersize',2,'color','w','linestyle','none');
  find(PSALqco=='3');  ho=line(X(ans),Y(ans)); set(ho,'marker','.','markersize',2,'color','w','linestyle','none');
  find(PSALqco=='4');  ho=line(X(ans),Y(ans)); set(ho,'marker','s','markersize',8,'color','w','linestyle','none');
  find(PSALqcn=='1');    hr=line(X(ans),Y(ans)); set(hr,'marker','d','markersize',8,'color','w','linestyle','none');
  find(PSALqcn=='4');    hn=line(X(ans),Y(ans)); set(hn,'marker','o','markersize',8,'color','w','linestyle','none');
  % find(PRESqco=='4'); ho=line(X(ans),Y(ans)); set(ho,'marker','s','markersize',8,'color',[.5 .5 .5],'linestyle','none');
  % find(PRESqcn=='1');   hr=line(X(ans),Y(ans)); set(hr,'marker','d','markersize',8,'color',[.5 .5 .5],'linestyle','none');
  % find(PRESqcn=='4');   hn=line(X(ans),Y(ans)); set(hn,'marker','o','markersize',8,'color',[.5 .5 .5],'linestyle','none');
  print(gcf,'-dpng',[outfiles{I}(1:end-4),'_HMsal.png']);
  figure(1003); axes(ap(2)); % Go to pressure plot
  % See above for RTQC flag markings on pressure plot.
  find(PRESqcn=='1');    hPn{1}=line(X(ans),Y(ans)); set(hPn{1},'marker','d','markersize',8,'color','r','linestyle','none','tag','qcflag','userdata','DMQC flag ''1''');
  find(PRESqcn=='4');    hPn{4}=line(X(ans),Y(ans)); set(hPn{4},'marker','o','markersize',8,'color','r','linestyle','none','tag','qcflag','userdata','DMQC flag ''4''');
  %
  try % Try, because sometimes there are no flags at all
    findobj(ap(2),'tag','qcflag'); legend(ans,get(ans,'userdata'),'location','best'); % Legend with any added flag marks (flexible legend using tag and userdata above)
  end
  % hfl=findobj(ap(2),'tag','qcflag')
  % legt=get(hfl,'userdata');
  % [hPo{2},hPo{3},hPo{4},hPo{8},hPo{9},hPn{1},hPn{4}]; % flag marks
  % if ~isempty(ans)
  %   legend(ans,'RTQC flag ''2''','RTQC flag ''3''','RTQC flag ''4''','RTQC flag ''8''','RTQC flag ''9''','DMQC flag ''1''','DMQC flag ''4''','location','best')
  % end
  print(gcf,'-dpng',[outfiles{I}(1:end-4),'_PRES.png']);  

  
  % %%%%%% Save for OWC and WRITE_D: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Add 360 to match the reference data:
  LONG<0; LONG(ans)=LONG(ans)+360; 
  %whos LAT LONG DATES PRES SAL TEMP PTMP PROFILE_NO CYCLE_NUMBER 
  % It is necessary to name the adjusted variables for OWC. These
  % will be renamed to their proper names with _ADJUSTED in WRITE_D:
  PRES=PRES_ADJUSTED; TEMP=TEMP_ADJUSTED; SAL=PSAL_ADJUSTED;
  
  % We send in the salinity profiles previously flagged as uncorrectable
  % to OWC in order to see the whole series in the plots from OWC. It
  % does make it easier to judge everything. We will have to do without
  % the CPcorr corrections, but compared to salinty drifts or offsets,
  % they are minute. Not sure if we should do OWC on cycles after PSAL
  % is bad and greylisted, but it is nice to show in the plots:
  if ASD<=n
    SAL(:,ASD:end)=SALo; 
    PTMP(:,ASD:end)=PTMPo; 
  end
  % For PSAL greylisted floats, PREPARE_FLOATS depends on finding the
  % calseries in set_calseries.m in order to keep the bad profiles in
  % the data to be used by OW_CALIBRATION.

  save(outfiles{I},...
       'LAT','LONG','DATES','PRES','SAL','TEMP','PTMP','PROFILE_NO',...
       'JULD','CYCLE_NUMBER',...
       'PROFILE_*_N','Rfiles','scientific_calib','PSAL_ADJUSTED_Cnew',...
       'ctdmodel','PSAL'); 
	% First row is for OWC, as well as OPERATOR_CPCOR_NEW and WRITE_D;
	% Second row is also used in write_D and ...?
	% Third row and PSAL is for WRITE_D;
	% Fourth row is for OPERATOR_CPCOR_NEW.
  
  % Save a local backup of your manual work, in float_working_dir since
  % outfiles are put in a non-backup/timemachine directory to avoid
  % performance issues and unnecessary use of backup space. It also
  % includes RTQC flags, since DMQC flags become RTQC flags in D-files.
  save(newflagfile,'*qco','*qcn'); disp(['prepare_floats: Writing RTQC and DMQC flags to ',newflagfile]);
  % For WRITE_D and future PREPARE_FLOATS.
  

  % [] Comments are excluded below, because it is hard to give comment during visual control.
  
  %%%%%%% Finalise comments and make the table of flags %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Go to that float's working dir and write snippets for qc-table:
  cd(float_working_dir);

  % Exclude the bad salinity-profiles from the stats:
  if ASD<=n
    PSALqco(:,ASD:end)=' '; PSALqcn(:,ASD:end)=' '; 
  end
  
  % if ~isempty(comments)
  %   comments(end-1:end)='. '; comments=['New qc=''4'' flags are due to: ',comments];  
  % else
  %   comments='No new qc=''4'' flags were found necessary.';
  % end
  snippet(comments);
  frm='%s & %s & %u & %u & %u & %s \\\\ \n'; % \multicolumn{2}{|c|}{Bottles}
  fid=fopen('flagtable.tex','w');
  fprintf(fid,'%s\n','\begin{tabular}{|ll|c|c|c|p{9cm}|}');
  fprintf(fid,'%s\n','\hline'); 
  %fprintf(fid,'%s \\\\ \n',['Variable & \multicolumn{5}{|c|}{RTQC}         & Reversed & New   & Affected cycles (cycle numbers)']);
  fprintf(fid,'%s \\\\ \n',['Variable & Flag & RTQC  & DMQC              & DMQC & Affected cycles (cycle numbers)']);
  fprintf(fid,'%s \\\\ \n',['         &      &       & reversed to ''1'' & new  &                                ']);
  fprintf(fid,'%s\n','\hline'); 
  fprintf(fid,'%s\n',' \hline');
  %                              sum of old flags        sum of reversed flags                sum of new flags        string with affected cycles
  fprintf(fid,frm,'POS ','''2''',sum(POSqco =='2','all'),sum(POSqco =='2'&POSqcn =='1','all'),sum(POSqcn =='2','all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(POSqco =='2'|POSqco =='2'&POSqcn =='1'|POSqcn =='2',1))))));
  fprintf(fid,frm,'    ','''3''',sum(POSqco =='3','all'),sum(POSqco =='3'&POSqcn =='1','all'),sum(POSqcn =='3','all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(POSqco =='3'|POSqco =='3'&POSqcn =='1'|POSqcn =='3',1))))));
  fprintf(fid,frm,'    ','''4''',sum(POSqco =='4','all'),sum(POSqco =='4'&POSqcn =='1','all'),sum(POSqcn =='4','all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(POSqco =='4'|POSqco =='4'&POSqcn =='1'|POSqcn =='4',1))))));
  fprintf(fid,frm,'    ','''8''',sum(POSqco =='8','all'),sum(POSqco =='8'&POSqcn =='1','all'),sum(POSqcn =='8','all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(POSqco =='8'|POSqco =='8'&POSqcn =='1'|POSqcn =='8',1))))));
  fprintf(fid,frm,'    ','''9''',sum(POSqco =='9','all'),sum(POSqco =='9'&POSqcn =='1','all'),sum(POSqcn =='9','all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(POSqco =='9'|POSqco =='9'&POSqcn =='1'|POSqcn =='9',1))))));
  fprintf(fid,'%s\n','\hline'); 
  fprintf(fid,frm,'JULD','''2''',sum(JULDqco=='2','all'),sum(JULDqco=='2'&JULDqcn=='1','all'),sum(JULDqcn=='2','all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(JULDqco=='2'|JULDqco=='2'&JULDqcn=='1'|JULDqcn=='2',1))))));
  fprintf(fid,frm,'    ','''3''',sum(JULDqco=='3','all'),sum(JULDqco=='3'&JULDqcn=='1','all'),sum(JULDqcn=='3','all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(JULDqco=='3'|JULDqco=='3'&JULDqcn=='1'|JULDqcn=='3',1))))));
  fprintf(fid,frm,'    ','''4''',sum(JULDqco=='4','all'),sum(JULDqco=='4'&JULDqcn=='1','all'),sum(JULDqcn=='4','all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(JULDqco=='4'|JULDqco=='4'&JULDqcn=='1'|JULDqcn=='4',1))))));
  fprintf(fid,frm,'    ','''8''',sum(JULDqco=='8','all'),sum(JULDqco=='8'&JULDqcn=='1','all'),sum(JULDqcn=='8','all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(JULDqco=='8'|JULDqco=='8'&JULDqcn=='1'|JULDqcn=='8',1))))));
  fprintf(fid,frm,'    ','''9''',sum(JULDqco=='9','all'),sum(JULDqco=='9'&JULDqcn=='1','all'),sum(JULDqcn=='9','all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(JULDqco=='9'|JULDqco=='9'&JULDqcn=='1'|JULDqcn=='9',1))))));
  fprintf(fid,'%s\n','\hline'); 
  fprintf(fid,frm,'PRES','''2''',sum(PRESqco=='2','all'),sum(PRESqco=='2'&PRESqcn=='1','all'),sum(PRESqcn=='2','all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(PRESqco=='2'|PRESqco=='2'&PRESqcn=='1'|PRESqcn=='2',1))))));
  fprintf(fid,frm,'    ','''3''',sum(PRESqco=='3','all'),sum(PRESqco=='3'&PRESqcn=='1','all'),sum(PRESqcn=='3','all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(PRESqco=='3'|PRESqco=='3'&PRESqcn=='1'|PRESqcn=='3',1))))));
  fprintf(fid,frm,'    ','''4''',sum(PRESqco=='4','all'),sum(PRESqco=='4'&PRESqcn=='1','all'),sum(PRESqcn=='4','all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(PRESqco=='4'|PRESqco=='4'&PRESqcn=='1'|PRESqcn=='4',1))))));
  fprintf(fid,frm,'    ','''9''',sum(PRESqco=='9','all'),sum(PRESqco=='9'&PRESqcn=='1','all'),sum(PRESqcn=='9','all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(PRESqco=='9'|PRESqco=='9'&PRESqcn=='1'|PRESqcn=='9',1))))));
  fprintf(fid,'%s\n','\hline'); 
  fprintf(fid,frm,'PSAL','''2''',sum(PSALqco=='2','all'),sum(PSALqco=='2'&PSALqcn=='1','all'),sum(PSALqcn=='2','all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(PSALqco=='2'|PSALqco=='2'&PSALqcn=='1'|PSALqcn=='2',1))))));
  fprintf(fid,frm,'    ','''3''',sum(PSALqco=='3','all'),sum(PSALqco=='3'&PSALqcn=='1','all'),sum(PSALqcn=='3','all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(PSALqco=='3'|PSALqco=='3'&PSALqcn=='1'|PSALqcn=='3',1))))));
  fprintf(fid,frm,'    ','''4''',sum(PSALqco=='4','all'),sum(PSALqco=='4'&PSALqcn=='1','all'),sum(PSALqcn=='4','all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(PSALqco=='4'|PSALqco=='4'&PSALqcn=='1'|PSALqcn=='4',1))))));
  fprintf(fid,frm,'    ','''9''',sum(PSALqco=='9','all'),sum(PSALqco=='9'&PSALqcn=='1','all'),sum(PSALqcn=='9','all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(PSALqco=='9'|PSALqco=='9'&PSALqcn=='1'|PSALqcn=='9',1))))));
  fprintf(fid,'%s\n','\hline'); 
  fprintf(fid,frm,'TEMP','''2''',sum(TEMPqco=='2','all'),sum(TEMPqco=='2'&TEMPqcn=='1','all'),sum(TEMPqcn=='2','all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(TEMPqco=='2'|TEMPqco=='2'&TEMPqcn=='1'|TEMPqcn=='2',1))))));
  fprintf(fid,frm,'    ','''3''',sum(TEMPqco=='3','all'),sum(TEMPqco=='3'&TEMPqcn=='1','all'),sum(TEMPqcn=='3','all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(TEMPqco=='3'|TEMPqco=='3'&TEMPqcn=='1'|TEMPqcn=='3',1))))));
  fprintf(fid,frm,'    ','''4''',sum(TEMPqco=='4','all'),sum(TEMPqco=='4'&TEMPqcn=='1','all'),sum(TEMPqcn=='4','all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(TEMPqco=='4'|TEMPqco=='4'&TEMPqcn=='1'|TEMPqcn=='4',1))))));
  fprintf(fid,frm,'    ','''9''',sum(TEMPqco=='9','all'),sum(TEMPqco=='9'&TEMPqcn=='1','all'),sum(TEMPqcn=='9','all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(TEMPqco=='9'|TEMPqco=='9'&TEMPqcn=='1'|TEMPqcn=='9',1))))));
  fprintf(fid,'%s\n','\hline'); 
  fprintf(fid,'%s\n','\end{tabular}');
  fclose(fid);
  
  
  
end	% Loop floats

cd(my_working_dir); % Go back to the main working dir

% For a final check of all profiles after this:
% figure
% for j=1:107; plot(PSAL(:,j),-PRES(:,j),'-',PSAL_ADJUSTED(:,j),-PRES(:,j),'--',SAL(:,j),-PRES(:,j),':'); input(['Profile ',int2str(j),'. Next?']); end



% Checklist on variables for OWC:
% LAT		√ (1×n, in decimal degrees, −ve means south of the equator, e.g. 20.5S = −20.5)
% LONG		√ (1×n, in decimal degrees, from 0 to 360, e.g. 98.5W in the eastern Pacific = 261.5E)
% DATES 	√ (1×n, in decimal year, e.g. 10 Dec 2000 = 2000.939726)
%		  *Note that this date format is different from that used in the reference database.
% PRES		√ (m×n, in dbar, monotonically increasing; i.e. 1st element of the column is the
%		  shallowest pressure, and subsequent values are unique and increasing)
% SAL		√ (m×n, in PSS-78)
% TEMP		√ (m×n, in-situ ITS-90)
% PTMP  	√ (m×n, potential temperature referenced to zero pressure)
% PROFILE_NO	√ (1×n, this can go from 0 to n−1, or 1 to n)
%		√ m = maximum number of observed levels from the float
%		√ n = number of profiles in the float time series




% DECIMATE:
% std_levels =[ 0          10          20          30          50 ...
% 	      75 100         125         150         200         250 ...
% 	      300 400         500         600         700         800 ...
% 	      900 1000        1100        1200        1300        1400 ...
% 	      1500 1750        2000        2200        2500 ...
% 	      3000     3500]';  % could just use standard levels for vertical grid
%floor(linspace(1,m,10)); i=setdiff(1:m,ans);
%[PRES,PTMP,SAL,TEMP]=delrow(i,PRES,PTMP,SAL,TEMP);

% 56:61  "No good climatological data found", but the wmo square dmqc
% files are there. 

% Some hacks: 
% On all figures in plot_diagnostics.m, frontalConstraintSAF.m, and anom.m:
% set(gcf,'OuterPosition',get(0,'ScreenSize')); % Even's addition

% Useful lines that can be inserted above for debugging purposes:
% disp('Post merge flags:')
% disp('JULDqco|JULDqcn; POSqco|POSqcn'), [JULDqco,'|',JULDqcn; POSqco,'|',POSqcn]
% disp('PRESqco|PRESqct|PRESqcn||PSALqco|PSALqct|PSALqcn'), repmat('|',m,1); [PRESqco,ans,PRESqct,ans,PRESqcn,ans,ans,PSALqco,ans,PSALqct,ans,PSALqcn]

% disp('From nc-file:')
% disp('JULDqco||JULDqco; POSqco||POSqco'), [JULDqco,'||',JULDqco; POSqco,'||',POSqco]
% disp('PRESqco||PSALqco'); repmat('|',m,1); [PRESqco,ans,ans,PSALqco]

% % Look at D-file QC, saved RTQC, and saved DMQC:
% disp('JULDqcob|JULDqco|JULDqcn; POSqcob|POSqco|POSqcn'), [JULDqcob,'|',JULDqco,'|',JULDqcn; POSqcob,'|',POSqco,'|',POSqcn]
% disp('PRESqcob|PRESqco|PRESqcn||PSALqcob|PSALqco|PSALqcn'), repmat('|',m,1); [PRESqcob,ans,PRESqco,ans,PRESqcn,ans,ans,PSALqcob,ans,PSALqco,ans,PSALqcn]

% disp('Post auto checks:')
% disp('JULDqco|JULDqcn; POSqco|POSqcn'), [JULDqco,'|',JULDqcn; POSqco,'|',POSqcn]
% disp('PRESqco|PRESqct|PRESqcn||PSALqco|PSALqct|PSALqcn'), repmat('|',m,1); [PRESqco,ans,PRESqct,ans,PRESqcn,ans,ans,PSALqco,ans,PSALqct,ans,PSALqcn]
    
% disp('Post replace ''X'':')
% disp('JULDqco|JULDqcn; POSqco|POSqcn'), [JULDqco,'|',JULDqcn; POSqco,'|',POSqcn]
% disp('PRESqco|PRESqct|PRESqcn||PSALqco|PSALqct|PSALqcn'), repmat('|',m,1); [PRESqco,ans,PRESqct,ans,PRESqcn,ans,ans,PSALqco,ans,PSALqct,ans,PSALqcn]
    
% disp('Post qcn asignment:')
% disp('JULDqco|JULDqcn; POSqco|POSqcn'), [JULDqco,'|',JULDqcn; POSqco,'|',POSqcn]
% disp('PRESqco|PRESqct|PRESqcn||PSALqco|PSALqct|PSALqcn'), repmat('|',m,1); [PRESqco,ans,PRESqct,ans,PRESqcn,ans,ans,PSALqco,ans,PSALqct,ans,PSALqcn]

% disp('Pre save:')
% disp('JULDqco|JULDqcn; POSqco|POSqcn'), [JULDqco,'|',JULDqcn; POSqco,'|',POSqcn]
% disp('PRESqco|PRESqcn||PSALqco|PSALqcn'), repmat('|',m,1); [PRESqco,ans,PRESqcn,ans,ans,PSALqco,ans,PSALqcn]
