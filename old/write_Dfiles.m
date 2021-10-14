% Writing D files. By Birgit Klein, J. Even Ø. Nilsen et al.
%
% This file is reading an exiting R\D-file created by Coriolis from
% rootdirin and writes D-files to be submitted to Coriolis to rootdirout.
% 
% It includes dmqc decisions for floats with no correction PSAL=PSAL_ADJUSTED,
% TEMP=TEMP_ADJUSTED, PRES=PRES_ADJUSTED
% or 
% with corrections PSAL=PSAL_ADJUSTED+ds, TEMP=TEMP_ADJUSTED, PRES=PRES_ADJUSTED.
% 
% The delayed mode is only performed on the primary profile and not on the
% NST profile, since according to the Argo User handbook it is allowed to
% perform the delayed mode at different times. According to the handbook
% the primary profiles should always be first, but check anyway.

init_dmqc;	% Load the setup. Float names, paths, etc. are set there. 
% float_names	is a cell object with the float numbers you have DMQC'ed.
% corr		is a corresponding cell with logical true for floats with
%		correction and false otherwise.
% Loading of float data as well as corresponding R-files is also done
% by init_dmqc (with on/off switch).

% One could imagine setting up a set of parameters for correction and
% reporting during owc analyses, and let this script run everything
% based on that, in a loop of I on all your floats. This is why the
% I is used throughout already now. Maybe later.	-Even -e
%
% To work on single float, simply change the loop statement:
%
for I=4				% Only one float
%for I=1:length(float_names)	% Loop all floats 

  % LOAD CALIBRATION RESULTS:
  load([float_dir,filesep,float_names{I}]); %load 6901815
  load([calib_dir,filesep,'cal_',float_names{I}]); % load cal_6901815
  ds=cal_SAL-SAL;			% The ds at all observations
  corr_ds=nanmean(ds);			% The mean ds for each cycle 
  corr_err=nanmean(cal_SAL_err);	% The mean error of for each cycle 
  lsnan=find(any(isnan(corr_ds)));	% Columns with NaNs for ds
  if ~isempty(lsnan)			% Check for NaNs
    disp(['Check these cycles for NaN corr term:',int2str(CYCLE_NUMBER(lsnan))]);
  end
  % do something about the missing corrections for the cycles listed -b
  % add code here, find out tomorrow -b

  b=dir([rootdirin{I},filesep,'*.nc']); % list all files in your rootout to work on all cycles 

  % [] it is exclusively profile files now, so no need to specify the 'D'. ??? 

  try, mkdir(rootdirout{I}); catch, end	% Make D-file directory
					
  for i=1:length(b) % Loop R-files
    % Read required information from exiting R-file created by Coriolis ion
    filenamein=[rootdirin{I},b(i).name];
    
    % filenamein(1)='R'; % [] Why this? 
  
    % WRITE_DFILES is reading an existing R-file created by Coriolis and includes
    % dmqc decisions from the OW method for floats which need no correction.
    
    % One R-file at a time is copied to a D-file:
    filenameout=b(i).name; filenameout(1)='D'; filenameout=[rootdirout{I},filenameout];
    system(['cp ',filenamein,' ',filenameout]);
    
    % [] Most likely all the rest could be done reading and writing
    % to the now copied outfile, without opening the infile. ???
    
    % INQUIRE SETTINGS FOR NLEVEL,NPROF,NCALIB, NHISTORY FROM THIS R-FILE:
    fidin=netcdf.open(filenamein,'NC_NOWRITE');	% Open the R-file
    pres=netcdf.getVar(fidin,netcdf.inqVarID(fidin,'PRES'));	% read pres matrix 
    [nlevel,nprof]=size(pres);					% Use it to determine dimensions
    if nprof>2
      display('nprof>2, please be aware of multiple profiles in this cycle, program only updating the primary profile in nprof=1')
    end
    % check if the primary sampling scheme is properly named in prof=1:
    vertical_sampling_scheme=netcdf.getVar(fidin,netcdf.inqVarID(fidin,'VERTICAL_SAMPLING_SCHEME')); 
    text=vertical_sampling_scheme(:,1)';
    wahr=strfind(text,'Primary');
    if ~wahr, display('Please check why no primary sampling!'); end
    % Determine the nparam dimension and ncalib status, here PARAMETER is used
    % which should have dimension STRING16,N_PARAM,N_CALIB,N_PROF
    parameter=netcdf.getVar(fidin,netcdf.inqVarID(fidin,'PARAMETER'));
    bb=size(parameter);
    nparam=bb(2);
    % ncalib should be stored in bb(3) in v3.1 files, but maybe not in older files:
    if length(bb)==2,	ncalib=1;
    else		ncalib=bb(3); end
    clear parameter bb a 
    % Find dimension of previous history field:
    history_start_pres=netcdf.getVar(fidin,netcdf.inqVarID(fidin,'HISTORY_START_PRES'));
    c=size(history_start_pres);
    nhistold=c(2);
    clear history_start_pres c
    % Determine time of creation of previous version:
    % [] Previous? This is now.
    dtext=datestr(now,30); datetext=[dtext(1:8) dtext(10:15)];
    % Get old time of creation for global atribute history:
    dateold=netcdf.getVar(fidin,netcdf.inqVarID(fidin,'DATE_CREATION'));

    
    % UPDATE OF THE D-FILE CAN NOW BEGIN:
    fid = netcdf.open(filenameout,'WRITE');

    %%%%%% HERE IS A LIST OF WHAT NEEDS TO BE UPDATED: %%%%%%%%%%%%%%%%%%%%%%
    % Birgit's list.
    % Update all variables and write them into new D-file.
    % The calculations and updating will be done in groups below this part.
    %
    % Here, √ means it is coded below; [] is an empty checkbox meaning it is pending.
    %
    % √ PRES_ADJUSTED -> PRES_ADJUSTED=PRES
    % √ PRES_ADJUSTED_QC -> PRES_ADJUSTED_QC=PRES_QC
    % PRES_ADJUSTED_ERROR -> 2.4 dbar
    % [] remember change data with bad qc in adjusted fields to fill values qc=4
    % [] remember to recalculate the PROFILE_<PARAM>_QC
    %
    % [] PROFILE_PRES_QC-> see argo dmqc manual, page 52 this in dmqc is
    % calculated from the adjusted_qc
    %
    % [] Please remember that APex floats need a surface pressure correction with
    % the correct cycle_numer -> PRES_ADJUSTED=PRES-surfacePRES 
    % has to be applied to PRES_ADJUSTED and then later on PSAL has to be
    % recalculated and to be put PSAL_ADJUSTED
    %
    % √ TEMP_ADJUSTED -> TEMP_ADJUSTED=TEMP;
    % √ TEMP_ADJUSTED_QC -> TEMP_ADJUSTED_QC-> TEMP_QC
    % √ TEMP_ADJUSTED_ERROR ->0.002
    % [] PROFILE_TEMP_QC -> recalculate from TEMP_ADJUSTED_QC
    % 
    % PSAL_ADJUSTED PSAL_ADJUSTED=PSAL -> then decide if correction is neeeded
    % corr=1 PSAL_ADJUSTED=PSAL+corr_ds
    % corr=0 PSAL_ADJUSTED=PSAL
    % [] But how to match up the cal-values from the float file and
    % the R-files? The number of cycles do not match.
    %
    % √ PSAL_ADJUSTED_QC-> PSAL_ADJUSTED_QC=PSAL_QC
    % [] PSAL_ADJUSTED_ERROR -> max of 0.01 and corr_err
    % PROFILE_PSAL_QC -> recalculate from PSAL_ADJUSTED_QC
    %
    %
    % remember change data with bad qc in adjusted fields to fill values
    % remember to recalculate the PROFILE_<PARAM>_QC
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    % PRES_ADJUSTED (35):
    % pres is loaded above
    netcdf.putVar(fid,netcdf.inqVarID(fid,'PRES_ADJUSTED'),pres);	
    
    % PRES_ADJUSTED_QC (36)
    netcdf.getVar(fid,netcdf.inqVarID(fid,'PRES_QC'));
    netcdf.putVar(fid,netcdf.inqVarID(fid,'PRES_ADJUSTED_QC'),ans);

    % PRES_ADJUSTED_ERROR (37): 
    % set pres_adjusted_error to 2.4 dbar, unless other information exist
    varid=netcdf.inqVarID(fid,'PRES_ADJUSTED_ERROR');		% var ID 
    pres_adjusted_error=netcdf.getVar(fid,varid);		% var values
    pres_adjusted_error==netcdf.getAtt(fid,varid,'_FillValue');	% Find fillvalues
    pres_adjusted_error(ans)=2.4;				% Set
    netcdf.putVar(fid,varid,pres_adjusted_error);		% Update

    % PROFILE_PRES_QC (28):
    % [] See argo dmqc manual, page 52 this in dmqc is calculated from the adjusted_qc
    % profile_pres_qc =
    % netcdf.putVar(fid,netcdf.inqVarID(fid,'PROFILE_PRES_QC'),profile_pres_qc);

    % TEMP_ADJUSTED (40):
    netcdf.getVar(fid,netcdf.inqVarID(fid,'TEMP'));
    netcdf.putVar(fid,netcdf.inqVarID(fid,'TEMP_ADJUSTED'),ans);
    
    % TEMP_ADJUSTED_QC (41):
    netcdf.getVar(fid,netcdf.inqVarID(fid,'TEMP_QC'));
    netcdf.putVar(fid,netcdf.inqVarID(fid,'TEMP_ADJUSTED_QC'),ans);
    
    % TEMP_ADJUSTED_ERROR (42):
    % set TEMP_ADJUSTED_ERROR to 0.002 °C, unless other information exist
    varid=netcdf.inqVarID(fid,'TEMP_ADJUSTED_ERROR');		% var ID 
    temp_adjusted_error=netcdf.getVar(fid,varid);		% var values
    temp_adjusted_error==netcdf.getAtt(fid,varid,'_FillValue');	% Find fillvalues
    temp_adjusted_error(ans)=0.002;				% Set
    netcdf.putVar(fid,varid,temp_adjusted_error);		% Update

    % PROFILE_TEMP_QC (29):
    % [] recalculate from TEMP_ADJUSTED_QC
    % profile_temp_qc=
    % netcdf.putVar(fid,netcdf.inqVarID(fid,'PROFILE_TEMP_QC'),profile_temp_qc);
    
    % PSAL_ADJUSTED (45)
    % [] There is a mismatch between R-files and cycles in the downloaded
    % netcdf.getVar(fid,netcdf.inqVarID(fid,'PSAL'));
    % float data. So there is a mismatch between cal
    % if corr{I},	ans(:,1)=ans(:,1)+corr_ds(i); end % Correct first column only
    % netcdf.putVar(fid,netcdf.inqVarID(fid,'PSAL_ADJUSTED'),ans);
    
    % PSAL_ADJUSTED_QC (46)
    netcdf.getVar(fid,netcdf.inqVarID(fid,'PSAL_QC'));
    netcdf.putVar(fid,netcdf.inqVarID(fid,'PSAL_ADJUSTED_QC'),ans);
  
    % PSAL_ADJUSTED_ERROR (47):
    % According to rules is either 0.01 or corr_err if that is greater than 0.01 
    % [] For Nordic Seas floats, we do not agree with this. 0.005 is our trust in the floats 
    % [] Also here comes the mismatch between R-file number (i) and cycle number.
    % [] Additional rule if correction corr_ds is >0.05 than either flag data as bad or increase error to min 0.015
    % varid=netcdf.inqVarID(fid,'PSAL_ADJUSTED_ERROR');		% var ID 
    % psal_adjusted_error=netcdf.getVar(fid,varid);		% var values
    % errfakt=max(0.010,abs(corr_err(i)))*sign(corr_err);	% Decide on error
    % psal_adjusted_error(:,1)=ones(nlevel,1)*errfakt;		% Adjust first column (only)
    % netcdf.putVar(fid,varid,psal_adjusted_error);		% Update

    % PROFILE_PSAL_QC (30)
    % [] How?
    % profile_ppsal_qc=
    % netcdf.putVar(fid,netcdf.inqVarID(fid,'PROFILE_PSAL_QC'),profile_psal_qc);
  
    % DATE_UPDATE (6):
    % Set to date when dmqc is performed, i.e., now (found above).
    netcdf.putVar(fid,netcdf.inqVarID(fid,'DATE_UPDATE'),datetext);

    % DATA_STATE_INDICATOR (15):
    % Follow rule in Argo dmqc manual (p.50) The variable DATA_STATE_INDICATOR
    % should record 2C or 2C+ and reference table 6 in user manual. 
    % No idea what 2C+ is?
    varid=netcdf.inqVarID(fidin,'DATA_STATE_INDICATOR');	% var ID 
    data_state_indicator=netcdf.getVar(fidin,varid);		% var values
    data_state_indicator(:,1)=['2C  ']';			% Adjust first column (only)
    netcdf.putVar(fid,varid,data_state_indicator);		% Update

    % DATA_MODE (16):
    varid=netcdf.inqVarID(fidin,'DATA_MODE');	% var ID 
    data_mode=netcdf.getVar(fidin,varid);	% var values
    data_mode(1)='D';				% Adjust first char (only)
    netcdf.putVar(fid,varid,data_mode);		% Update

    
    % SCIENTIFIC_CALIB_**** fields for floats without corrections
    % calib equation:
    pres_scientific_eq_dmqc='PRES_ADJUSTED = PRES';
    temp_scientific_eq_dmqc='TEMP_ADJUSTED = TEMP';
    psal_scientific_eq_dmqc='PSAL_ADJUSTED = PSAL';
    % calib_comment:
    pres_scientific_com_dmqc='No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar';
    psal_scientific_com_dmqc='No significant salinity drift detected. OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.';
    temp_scientific_com_dmqc='No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90';
    % coefficient:
    pres_scientific_coef_dmqc='none';
    temp_scientific_coef_dmqc='none';
    psal_scientific_coef_dmqc='none';
    if corr{I} % Change the fields necessary
      % SCIENTIFIC_CALIB_**** fields for floats with salinity corrections
      psal_scientific_eq_dmqc='PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r';
      psal_scientific_com_dmqc='Significant salinity drift present  - correction applied using OWC method (weighted least squares piecewise-fit).The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.';
      psal_scientific_coef_dmqc=['r= ' num2str(pcond_factor(i)) ', vertically averaged dS= ' num2str(corr_ds(i))];
    end

    % COPY THE NEW TEXTS INTO THE SCIENTIFIC FIELDS
    % 
    % First find sequence of parameters from field PARAMETER and find
    % index of primary profile (in case of NPROF>1) from field
    % VERTICAL_SAMPLING_SCHEME. 
    parameter=netcdf.getVar(fidin,netcdf.inqVarID(fidin,'PARAMETER')); 
    % Reihenfolger von T,S,P in den SCIENTIFIC Felder festlegen:
    % Assumes that both nprof have the same order because it is search on nprof=1
    text=cellstr(parameter(1:4,:,1,1)');
    it=find(strcmp(text,'TEMP')); if isempty(it), warning('no such parameter as TEMP, please check'); end
    is=find(strcmp(text,'PSAL')); if isempty(is), warning('no such parameter as PSAL, please check'); end
    ip=find(strcmp(text,'PRES')); if isempty(ip), warning('no such parameter as PRES, please check'); end
    % Bei Nprof>1, den Index für das primary sampling finden:
    vertical_sampling_scheme=netcdf.getVar(fidin,netcdf.inqVarID(fidin,'VERTICAL_SAMPLING_SCHEME')); 
    text=cellstr(vertical_sampling_scheme');
    in=find(contains(text,'Primary sampling:')); if isempty(in), warning('no entry for Primary sampling, please check'); end
    
    % Then update the fields in the D-files:
    
    % SCIENTIFIC_CALIB_EQUATION (61):
    % Für Nprof=1, rest remains emptyn (i.e., the first profile only)
    scientific_calib_equation=netcdf.getVar(fidin,netcdf.inqVarID(fidin,'SCIENTIFIC_CALIB_EQUATION')); 
    text256=repmat(' ',256,1);	% TEMP:
    text256(1:length(temp_scientific_eq_dmqc))=temp_scientific_eq_dmqc;
    scientific_calib_equation(:,it,ncalib,in)=text256;
    text256=repmat(' ',256,1);	% PSAL:
    text256(1:length(psal_scientific_eq_dmqc))=psal_scientific_eq_dmqc;
    scientific_calib_equation(:,is,ncalib,in)=text256;
    text256=repmat(' ',256,1);	% PRES:
    text256(1:length(pres_scientific_eq_dmqc))=pres_scientific_eq_dmqc;
    scientific_calib_equation(:,ip,ncalib,in)=text256;
    netcdf.putVar(fid,netcdf.inqVarID(fid,'SCIENTIFIC_CALIB_EQUATION'),scientific_calib_equation);
    clear scientific_calib_equation

    % SCIENTIFIC_CALIB_COEFFICIENT (62):
    scientific_calib_coefficient=netcdf.getVar(fidin,netcdf.inqVarID(fidin,'SCIENTIFIC_CALIB_COEFFICIENT'));
    text256=repmat(' ',256,1); % TEMP
    text256(1:length(temp_scientific_coef_dmqc))=temp_scientific_coef_dmqc;
    scientific_calib_coefficient(:,it,ncalib,in)=text256;
    text256=repmat(' ',256,1); % PSAL
    text256(1:length(psal_scientific_coef_dmqc))=psal_scientific_coef_dmqc;
    scientific_calib_coefficient(:,is,ncalib,in)=text256;
    text256=repmat(' ',256,1); % PRES
    text256(1:length(pres_scientific_coef_dmqc))=pres_scientific_coef_dmqc;
    scientific_calib_coefficient(:,ip,ncalib,in)=text256;
    netcdf.putVar(fid,netcdf.inqVarID(fid,'SCIENTIFIC_CALIB_COEFFICIENT'),scientific_calib_coefficient);
    clear scientific_calib_coefficient

    % SCIENTIFIC_CALIB_COMMENT (63):
    scientific_calib_comment=netcdf.getVar(fidin,netcdf.inqVarID(fidin,'SCIENTIFIC_CALIB_COMMENT')); 
    text256=repmat(' ',256,1);
    text256(1:length(temp_scientific_com_dmqc))=temp_scientific_com_dmqc;
    scientific_calib_comment(:,it,ncalib,in)=text256;
    text256=repmat(' ',256,1);
    text256(1:length(psal_scientific_com_dmqc))=psal_scientific_com_dmqc;
    scientific_calib_comment(:,is,ncalib,in)=text256;
    text256=repmat(' ',256,1);
    text256(1:length(pres_scientific_com_dmqc))=pres_scientific_com_dmqc;
    scientific_calib_comment(:,ip,ncalib,in)=text256;
    netcdf.putVar(fid,netcdf.inqVarID(fid,'SCIENTIFIC_CALIB_COMMENT'),scientific_calib_comment);
    clear scientific_calib_comment
    
    % SCIENTIFIC_CALIB_DATE (64):
    % date and dtext is found above
    scientific_calib_date=netcdf.getVar(fidin,netcdf.inqVarID(fidin,'SCIENTIFIC_CALIB_DATE')); 
    scientific_calib_date(:,it,ncalib,in)=datetext;
    scientific_calib_date(:,is,ncalib,in)=datetext;
    scientific_calib_date(:,ip,ncalib,in)=datetext;
    netcdf.putVar(fid,netcdf.inqVarID(fid,'SCIENTIFIC_CALIB_DATE'),scientific_calib_date);
    clear scientific_calib_date

    % NOW PREPARE TEXT FOR THE HISTORY FIELDS:
    % Text fields to be used according to reference tables 
    text_histinstitution='GE  '; % see reference table 4 in user manual
    text_histaction='IP  ';      % see reference table 7 in user manual
    text_histstep='ARSQ';        % see reference table 12 in user manual
    text_histsoftware='OWC ';    % locally defined
    text_histsoftwarerelease='2.0 '; % locally defined
    text_histreference='ARGO CTD ref. database: CTD_for_DMQC_2019V01 + ARGO_for_DMQC_2020V01';  % please update whenever a new version of ref data is available
    text_histdate=datetext;
    text_histparameter='PSAL            '; % see note on page 103 in user manual, should actually be fill value
    
    % [] Neither presanf nor presend is defined
    
    histstartpres=presanf;                 % should actually be fill value
    histendpres=presend;                   % should actually be fill value
    text_histqctest='                ';    % is fill value
    histpreviousvalue=99999;               % is fill value
    
    % HISTORY_INSTITUTION (48):
    % Es scheint das bei der ersten belegung eines unlimited feldes
    % start und count Felder braucht ansonsten bekommt man eine
    % fehlermeldung. 
    % The number of input elements does not match the variable size.
    history_institution=netcdf.getVar(fidin,netcdf.inqVarID(fidin,'HISTORY_INSTITUTION'));
    nhistnew=nhistold+1;
    history_institution_new=history_institution;
    history_institution_new(:,in,nhistnew)=text_histinstitution;
    netcdf.putVar(fid,netcdf.inqVarID(fid,'HISTORY_INSTITUTION'),[0,0,0],[4,nprof,nhistnew],history_institution_new);

    % HISTORY_STEP (49):
    history_step=netcdf.getVar(fidin,netcdf.inqVarID(fidin,'HISTORY_STEP'));
    history_step_new=history_step;
    history_step_new(:,in,nhistnew)=text_histstep;
    netcdf.putVar(fid,netcdf.inqVarID(fid,'HISTORY_STEP'),history_step_new);

    % HISTORY_SOFTWARE (50):
    history_software=netcdf.getVar(fidin,netcdf.inqVarID(fidin,'HISTORY_SOFTWARE'));
    history_software_new=history_software;
    history_software_new(:,in,nhistnew)=text_histsoftware;
    netcdf.putVar(fid,netcdf.inqVarID(fid,'HISTORY_SOFTWARE'),history_software_new);

    % HISTORY_SOFTWARE_RELEASE (51):
    history_software_release=netcdf.getVar(fidin,netcdf.inqVarID(fidin,'HISTORY_SOFTWARE_RELEASE'));
    history_software_release_new=history_software_release;
    history_software_release_new(:,in,nhistnew)=text_histsoftwarerelease;
    netcdf.putVar(fid,netcdf.inqVarID(fid,'HISTORY_SOFTWARE_RELEASE'),history_software_release_new);

    % HISTORY_REFERENCE (52):
    history_reference=netcdf.getVar(fidin,netcdf.inqVarID(fidin,'HISTORY_REFERENCE'));
    history_reference_new=history_reference;
    history_reference_new(:,in,nhistnew)=text_histreference;
    netcdf.putVar(fid,netcdf.inqVarID(fid,'HISTORY_REFERENCE'),history_reference_new);

    % HISTORY_DATE (53):
    history_date=netcdf.getVar(fidin,netcdf.inqVarID(fidin,'HISTORY_DATE'));
    history_date_new=history_date;
    history_date_new(:,in,nhistnew)=text_histdate;
    netcdf.putVar(fid,netcdf.inqVarID(fid,'HISTORY_DATE'),history_date_new);

    % HISTORY_ACTION (54):
    history_action=netcdf.getVar(fidin,netcdf.inqVarID(fidin,'HISTORY_ACTION'));
    history_action_new=history_action;
    history_action_new(:,in,nhistnew)=text_histaction;
    netcdf.putVar(fid,netcdf.inqVarID(fid,'HISTORY_ACTION'),history_action_new);

    % HISTORY_PARAMETER (55):
    history_parameter=netcdf.getVar(fidin,netcdf.inqVarID(fidin,'HISTORY_PARAMETER'));
    history_parameter_new=history_parameter;
    history_parameter_new(:,in,nhistnew)=text_histparameter;
    netcdf.putVar(fid,netcdf.inqVarID(fid,'HISTORY_PARAMETER'),history_parameter_new);

    % HISTORY_START_PRES (56):
    history_start_pres=netcdf.getVar(fidin,netcdf.inqVarID(fidin,'HISTORY_START_PRES'));
    history_start_pres_new=history_start_pres;
    history_start_pres_new(in,nhistnew)=presanf;
    if nprof>1, history_start_pres_new(2,nhistnew)=99999; end
    netcdf.putVar(fid,netcdf.inqVarID(fid,'HISTORY_START_PRES'),history_start_pres_new);

    % HISTORY_STOP_PRES (57):
    history_stop_pres=netcdf.getVar(fidin,netcdf.inqVarID(fidin,'HISTORY_STOP_PRES'));
    history_stop_pres_new=history_stop_pres;
    history_stop_pres_new(in,nhistnew)=presend;
    if nprof>1, history_stop_pres_new(2,nhistnew)=99999; end
    netcdf.putVar(fid,netcdf.inqVarID(fid,'HISTORY_STOP_PRES'),history_stop_pres_new);

    % HISTORY_PREVIOUS_VALUS (58):
    history_previous_value=netcdf.getVar(fidin,netcdf.inqVarID(fidin,'HISTORY_PREVIOUS_VALUE'));
    history_previous_value_new=history_previous_value;
    history_previous_value_new(in,nhistnew)=histpreviousvalue;
    if nprof>1, history_previous_value_new(2,nhistnew)=99999; end
    netcdf.putVar(fid,netcdf.inqVarID(fid,'HISTORY_PREVIOUS_VALUE'),history_previous_value_new);

    % HISTORY_QCTEST (59):
    history_qctest=netcdf.getVar(fidin,netcdf.inqVarID(fidin,'HISTORY_QCTEST'));
    history_qctest_new=history_qctest;
    history_qctest_new(:,in,nhistnew)=text_histqctest;
    netcdf.putVar(fid,netcdf.inqVarID(fid,'HISTORY_QCTEST'),history_qctest_new);

    % % finished, now close files
    % 
    % [] This final part is not yet checked or understood by me. -e
    %
    % reenter define mode
    % mochte den Aufruf nicht mit NC_GLOBAL musste das in nc_global umändern, 
    % Cannot find an exact (case-sensitive) match for 'NC_GLOBAL'
    %
    %The closest match is: nc_global in C:\Program
    %Files\MATLAB\R2013a\mexnc-2.0.30\nc_global.m
    %
    NC_GLOBAL = netcdf.getConstant('NC_GLOBAL');
    netcdf.reDef(fid)
    dateold=dateold';
    oldtext=[dateold(1:4) '-' dateold(5:6) '-' dateold(7:8) 'T' dateold(9:10) ':' dateold(11:12) ':' dateold(13:14) 'Z creation; '];
    updatetext=[oldtext datetext(1:4) '-' datetext(5:6) '-' datetext(7:8) 'T' datetext(9:10) ':' datetext(11:12) ':' datetext(13:14) 'Z last update (BSH ARSQ software)'];
    netcdf.putAtt(fid,NC_GLOBAL,'history',updatetext)
    netcdf.endDef(fid)
    netcdf.close(fid)
    netcdf.close(fidin)
    clear nlevel nprof ncalib nparam nhistold nhistnew
    clear dtext datetext dateold
  end
end
