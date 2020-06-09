% 
%
% This file is reading an exiting R\D-file created by Coriolis from
% rootdirin and writes D-files to be submitted to Coriolis to rootdirout
% 
% it includes dmqc decisions for floats with no correction PSAL=PSAL_ADJUSTED,
% TEMP=TEMP_ADJUSTED, PRES=PRES_ADJUSTED
% or 
% with corrections PSAL=PSAL_ADJUSTED+ds,
% TEMP=TEMP_ADJUSTED, PRES=PRES_ADJUSTED
%
% 
% the delayed mode is only performed on the primary profile and not on the
% NST profile, since according to the Argo User handbook it is allowed to
% perform the delayed mode at different times. According to the handbook
% the primary profiles should always be first, but check anyway.

init_dmqc;	% Float names, paths, etc. are set there. 


if logical(0) % ----- SWITCH ON IN ORDER TO DOWNLOAD ALL THE R-FILES FIRST ------------------------
  for I=5 %1:length(float_names)
    mkdir(rootdirin{I}); 
    system(['lftp -e ''lcd ',rootdirin{I},' ; mget R',float_names{I},'_*.nc ; bye'' ftp.ifremer.fr/ifremer/argo/dac/coriolis/',float_names{I},'/profiles/']);
    %mkdir(rootdirout{I}); 
  end
end % ----------------------------------------------------------------------------------------

% Treat one float at a time, but base everything on your own set of float names, as set in INIT_DMQC:
%float_names={'6903548' '6903549' '6903550' '6903551' '6903553' '6903554' '6903555' '6903557' '6903558' '6903559' '6903560' '6903561' '6903562' '6903563' '6903564'};
%
% One could imagine setting up a set of parameters for correction and
% reporting during owc analyses, and let this script run everything
% based on that, in an loop of I on all your floats. This is why the
% I is used throughout already now.		-Even -e
%
% select float
%
I=5; % Float number 5 is float 6903553.

%%float=float_names{I}; %6901815;
%%rootdirout='h:\DMQC_offline\DMQC_files\';
%%pathout=[rootdirout{I},num2str(float),filesep,'profiles',filesep];
%%cd (pathout) % No, stay in working directory and use full paths. -e
%
% determine if to correct corr=1 or not corr=0
% will still need a cycle number where to start with corrections -> use
% corr_ds to set it 0 for the part of the
% 
corr=0;

% needs output from your owc run the matfile created to perform the OWC
% run and the cal files that resulted from it, copy it in the same
% directory as the D-files to be updated. 
% No copying. Just load with paths. -e
%
load([float_dir,filesep,float_names{I}]); %load 6901815
load([calib_dir,filesep,'cal_',float_names{I}]); % load cal_6901815
ds=cal_SAL-SAL;
%
% calculate the mean ds for each cycle and the mean error of the mapping
% for each cycle
%
a=size(ds);
for i=1:a(2)
  l1=find(~isnan(cal_SAL(:,i)));
  corr_ds(i)=mean(ds(l1,i));
  corr_err(i)=mean(cal_SAL_err(l1,i));
end
lsnan=find(isnan(corr_ds));
if ~isempty(lsnan)
  disp('check cycles for NaN corr term')
  CYCLE_NUMBER(lsnan)
end
%
% do something about the missing corrections for the cycles listed
%
%
% add code here, find out tomorrow
%
%
%%b=dir('D*.nc'); % list all the d-files in your rootout to work on all cycles 

system(['cp -r ',rootdirin{I},' ',rootdirout{I}]); % Copy all files
b=dir([rootdirout{I},filesep,'*.nc']); % list all files in your rootout to work on all cycles 
% it is exclusively profile files now, so no need to specify the 'D'.

return
%
%
for i=1:length(b)
  filename=b(i).name;
  filenameout=filename;
  %
  % read required information from exiting R-file created by Coriolis ion
  %
  filenamein=filename;
  filenamein(1)='R';
  % This file is reading an exiting R-file created by Coriolis and includes
  % dmqc decisions from the OW method for floats which need no correction
  %
rootdirin='h:\DMQC_offline\';
pathin=[rootdirin num2str(float_names{I}) '\profiles\'];
fidinfile=netcdf.open(strcat(pathin,filenamein),'NC_NOWRITE');
%
% Inquire settings for nlevel,nprof,ncalib, nhistory from this file
%
varidin=netcdf.inqVarID(fidinfile,'PRES');
pres=netcdf.getVar(fidinfile,varidin);
a=size(pres);
nlevel=a(1);
nprof=a(2);
if nprof>2
    display('nprof>2, please be aware of multiple profiles in this cycle, program only updating the primary profile in nprof=1')
end
%
% Basic format checking in the existing file
%
%
% check if the primary sampling schem is properly named in nprof=1
%
varidin=netcdf.inqVarID(fidinfile,'VERTICAL_SAMPLING_SCHEME');
vertical_sampling_scheme=netcdf.getVar(fidinfile,varidin); 
text=vertical_sampling_scheme(:,1);
text=text';
wahr=strfind(text,'Primary');
if wahr~=1
    display('please check why no primary sampling')
end
%
% determine the nparam dimension and ncalib status, here PARAMETER is used
% which should have dimension STRING16,N_PARAM,N_CALIB,N_PROF
%
varidin=netcdf.inqVarID(fidinfile,'PARAMETER');
parameter=netcdf.getVar(fidinfile,varidin);
bb=size(parameter);
%
% 
%
nparam=bb(2);
%
% ncalib should be stored in bb(3) in v3.1 files, but maybe not in older
% files
%
if length(bb)==2
    ncalib=1;
else
    ncalib=bb(3);
end
clear parameter bb a 
%
% find dimension of previous history field
%
varidin=netcdf.inqVarID(fidinfile,'HISTORY_START_PRES');
history_start_pres=netcdf.getVar(fidinfile,varidin);
c=size(history_start_pres);
nhistold=c(2);
clear history_start_pres c
%
% determine time of creation of previous verion
%
dtext=datestr(now,30);
datetext=[dtext(1:8) dtext(10:15)];
%
% get old time of creation for global Atribute history
%
varidin=netcdf.inqVarID(fidinfile,'DATE_CREATION');
dateold=netcdf.getVar(fidinfile,varidin);
%
% Update of the D-file can now beginn
%
fid = netcdf.open(strcat(path1out,b(i).name),'WRITE');
%
% update all variables and write them into new D-file
%
% which fields need to be updated?
%
% PRES_ADJUSTED -> PRES_ADJUSTED=PRES
%  PRES_ADJUSTED_QC -> PRES_ADJUSTED_QC=PRES_QC
% PRES_ADJUSTED_ERROR -> 2.4 dbar
% remember change data with bad qc in adjusted fields to fill values qc=4
% remember to recalculate the PROFILE_<PARAM>_QC
%
% PROFILE_PRES_QC-> see argo dmqc manual, page 52 this in dmqc is
% calculated from the adjusted_qc
%
% Please remember that APex floats need a surface pressure correction with
% the correct cycle_numer -> PRES_ADJUSTED=PRES-surfacePRES 
% has to be applied to PRES_ADJUSTED and then later on PSAL has to be
% recalculated and to be put PSAL_ADJUSTED
%
%
% TEMP_ADJUSTED -> TEMP_ADJUSTED=TEMP;
% TEMP_ADJUSTED_QC -> TEMP_ADJUSTED_QC-> TEMP_QC
% TEMP_ADJUSTED_ERROR ->0.002
% PROFILE_TEMP_QC -> recalculate from TEMP_ADJUSTED_QC
% 
% PSAL_ADJUSTED PSAL_ADJUSTED=PSAL -> then decide if correction is neeeded
% corr=1 PSAL_ADJUSTED=PSAL+corr_ds
% corr=0 PSAL_ADJUSTED=PSAL
% PSAL_ADJUSTED_QC-> PSAL_ADJUSTED_QC=PSAL_QC
% PSAL_ADJUSTED_ERROR -> max of 0.01 and corr_err
% PROFILE_PSAL_QC -> recalculate from PSAL_ADJUSTED_QC
%
%
% remember change data with bad qc in adjusted fields to fill values
% remember to recalculate the PROFILE_<PARAM>_QC
%
%
% PRES_ADJUSTED (35)
%
varidout=netcdf.inqVarID(fid,'PRES_ADJUSTED');
netcdf.putVar(fid,varidout,pres_adjusted);
%
% PRES_ADJUSTED_QC (36)
%
varidout=netcdf.inqVarID(fid,'PRES_ADJUSTED_QC');
netcdf.putVar(fid,varidout,pres_adjusted_qc);
%
% PRES_ADJUSTED_ERROR (37)
%
% set PRES_ADJUSTED_ERROR to 2.4 dbar, unless other information exist
%
varidout=netcdf.inqVarID(fid,'PRES_ADJUSTED_ERROR');
netcdf.putVar(fid,varidout,pres_adjusted_error);
%
% PROFILE_PRES_QC (28)
%
varidout=netcdf.inqVarID(fid,'PROFILE_PRES_QC');
netcdf.putVar(fid,varidout,profile_pres_qc);
%
% TEMP_ADJUSTED (40)
%
varidout=netcdf.inqVarID(fid,'TEMP_ADJUSTED');
netcdf.putVar(fid,varidout,temp_adjusted);
%
% TEMP_ADJUSTED_QC (41)
%
varidout=netcdf.inqVarID(fid,'TEMP_ADJUSTED_QC');
netcdf.putVar(fid,varidout,temp_adjusted_qc);
%
% TEMP_ADJUSTED_ERROR (42)
% set TEMP_ADJUSTED_ERROR to 0.002 °C, unless other information exist
%
varidout=netcdf.inqVarID(fid,'TEMP_ADJUSTED_ERROR');
netcdf.putVar(fid,varidout,temp_adjusted_error);
%
% PROFILE_TEMP_QC (29)
%
varidout=netcdf.inqVarID(fid,'PROFILE_TEMP_QC');
netcdf.putVar(fid,varidout,profile_temp_qc);
%
%
% PSAL_ADJUSTED (45)
%
varidout=netcdf.inqVarID(fid,'PSAL_ADJUSTED');
netcdf.putVar(fid,varidout,psal_adjusted);
%
% PROFILE_PSAL_QC (30)
%
varidout=netcdf.inqVarID(fid,'PROFILE_PSAL_QC');
netcdf.putVar(fid,varidout,profile_psal_qc);
%
% PSAL_ADJUSTED_QC (46)
%
varidout=netcdf.inqVarID(fid,'PSAL_ADJUSTED_QC');
netcdf.putVar(fid,varidout,psal_adjusted_qc);
%
% PSAL_ADJUSTED_ERROR (47) 
% according to rules is either 0.01 or corr_err if that is greater then
% 0.01 
% Additional rule if correction corr_ds is >0.05 than either flag data as
% bad or increase error to min 0.015
%
%
errfakt=0.010;
%
% ersetze Standardfehler durch mapping fehler, wenn >0.01
%
if corr_err(i)>0.01
    errfakt=corr_err(i);
end
psal_adjusted_error(:,1)=ones(nlevel,1)*errfakt;
varidout=netcdf.inqVarID(fid,'PSAL_ADJUSTED_ERROR');
netcdf.putVar(fid,varidout,psal_adjusted_error);
%
% DATE_UPDATE (6)
% set this to date when dmqc is performed
%
date_update=datetext;
varidout=netcdf.inqVarID(fid,'DATE_UPDATE');
netcdf.putVar(fid,varidout,date_update);
clear date_update varidin
%
%
% DATA_STATE_INDICATOR (15)
%
% varidin=netcdf.inqVarID(fidinfile,'DATA_STATE_INDICATOR');
varidin=netcdf.inqVarID(fidinfile,'DATA_STATE_INDICATOR');
data_state_indicator=netcdf.getVar(fidinfile,varidin);
%
% follow rule in Argo dmqc manual (p.50) The variable DATA_STATE_INDICATOR
% should record ‘2C’ or ‘2C+’ and reference table 6 in user manual. No idea
% what 2C+ is?
%
data_state_indicator(:,1)='2C  ';
varidout=netcdf.inqVarID(fid,'DATA_STATE_INDICATOR');
netcdf.putVar(fid,varidout,data_state_indicator);
%
% DATA_MODE (16)
%
% varidin=netcdf.inqVarID(fidinfile,'DATA_MODE');
varidin=netcdf.inqVarID(fidinfile,'DATA_MODE');
data_mode=netcdf.getVar(fidinfile,varidin);
data_mode(1)='D';
varidout=netcdf.inqVarID(fid,'DATA_MODE');
netcdf.putVar(fid,varidout,data_mode);
if corr==1
%
% SCIENTIFIC_CALIB_**** fields for floats with corrections
% 
% 
% calib equation
%
pres_scientific_eq_dmqc='PRES_ADJUSTED = PRES';
temp_scientific_eq_dmqc='TEMP_ADJUSTED = TEMP';
psal_scientific_eq_dmqc='PSAL + dS, where dS is calculated from a potential conductivity (ref to 0 dbar) multiplicative adjustment term r';
%
% calib_comment
%
pres_scientific_com_dmqc='No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar';
psal_scientific_com_dmqc='Significant salinity drift present  - correction applied using OWC method (weighted least squares piecewise-fit).The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.';
temp_scientific_com_dmqc='No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90';
%
% coefficient
%
pres_scientific_coef_dmqc='none';
temp_scientific_coef_dmqc='none';
psal_scientific_coef_dmqc=['r= ' num2str(pcond_factor(i)) ', vertically averaged dS= ' num2str(corr_ds(i))];
end
if corr==0
%
% SCIENTIFIC_CALIB_**** fields for floats without corrections
% 
% 
% calib equation
%
pres_scientific_eq_dmqc='PRES_ADJUSTED = PRES';
temp_scientific_eq_dmqc='TEMP_ADJUSTED = TEMP';
psal_scientific_eq_dmqc='PSAL_ADJUSTED = PSAL';
%
% calib_comment
%
pres_scientific_com_dmqc='No significant pressure drift detected. Calibration error is manufacturer specified accuracy in dbar';
psal_scientific_com_dmqc='No significant salinity drift detected. OW method (weighted least squares fit) adopted. The quoted error is max[0.01, 1xOW uncertainty] in PSS-78.';
temp_scientific_com_dmqc='No significant temperature drift detected. Calibration error is manufacturer specified accuracy with respect to ITS-90';
%
% coefficient
%
pres_scientific_coef_dmqc='none';
temp_scientific_coef_dmqc='none';
psal_scientific_coef_dmqc='none';
end
%
% copy the new texts into the SCIENTIFIC fields
% 
% first find sequence of parameters from field PARAMETER
% and find index of primary profile (in case of NPROF>1)
% from field VERTICAL_SAMPLING_SCHEME
% 
varidin=netcdf.inqVarID(fidinfile,'PARAMETER');
parameter=netcdf.getVar(fidinfile,varidin); 
varidin=netcdf.inqVarID(fidinfile,'VERTICAL_SAMPLING_SCHEME');
vertical_sampling_scheme=netcdf.getVar(fidinfile,varidin); 
%
% Reihenfolger von T,S,P in den SCIENTIFIC Felder festlegen
% Assumes that both nprof have the same order because it is search on
% nprof=1
%
a=size(parameter);
for i=1:a(2)
text=parameter(:,i,1,1);
text=text';
text1=text(1:4);
tft(i)=strcmp('TEMP',text1);
tfs(i)=strcmp('PSAL',text1);
tfp(i)=strcmp('PRES',text1);
end
it=find(tft==1);
if isempty(it)
    warning('no such parameter as TEMP, please check')
end
is=find(tfs==1);
if isempty(is)
    warning('no such parameter as PSAL, please check')
end
ip=find(tfp==1);
if isempty(ip)
    warning('no such parameter as PRES, please check')
end
%
% Bei Nprof>1, den Index für das primary sampling finden
%
a=size(vertical_sampling_scheme);
for i=1:a(2)
text=vertical_sampling_scheme(:,i);
text=text';
text1=text(1:17);
tfn(i)=strcmp('Primary sampling:',text1);   
end
in=find(tfn==1);
if isempty(in)
    warning('no entry for Primary sampling, please check')
end
%
% SCIENTIFIC_CALIB_EQUATION (61), für Nprof=1, rest remains emptyn 
%
varidin=netcdf.inqVarID(fidinfile,'SCIENTIFIC_CALIB_EQUATION');
scientific_calib_equation=netcdf.getVar(fidinfile,varidin); 
text256=repmat(' ',256,1);
text256(1:length(temp_scientific_eq_dmqc))=temp_scientific_eq_dmqc;
scientific_calib_equation(:,it,ncalib,in)=text256;
%
text256=repmat(' ',256,1);
text256(1:length(psal_scientific_eq_dmqc))=psal_scientific_eq_dmqc;
scientific_calib_equation(:,is,ncalib,in)=text256;
%
text256=repmat(' ',256,1);
text256(1:length(pres_scientific_eq_dmqc))=pres_scientific_eq_dmqc;
scientific_calib_equation(:,ip,ncalib,in)=text256;
varidout=netcdf.inqVarID(fid,'SCIENTIFIC_CALIB_EQUATION');
netcdf.putVar(fid,varidout,scientific_calib_equation);
clear scientific_calib_equation varidin
%
% SCIENTIFIC_CALIB_COEFFICIENT (62)
%
varidin=netcdf.inqVarID(fidinfile,'SCIENTIFIC_CALIB_COEFFICIENT');
scientific_calib_coefficient=netcdf.getVar(fidinfile,varidin); 
text256=repmat(' ',256,1);
text256(1:length(temp_scientific_coef_dmqc))=temp_scientific_coef_dmqc;
scientific_calib_coefficient(:,it,ncalib,in)=text256;
%
text256=repmat(' ',256,1);
text256(1:length(psal_scientific_coef_dmqc))=psal_scientific_coef_dmqc;
scientific_calib_coefficient(:,is,ncalib,in)=text256;
%
text256=repmat(' ',256,1);
text256(1:length(pres_scientific_coef_dmqc))=pres_scientific_coef_dmqc;
scientific_calib_coefficient(:,ip,ncalib,in)=text256;
varidout=netcdf.inqVarID(fid,'SCIENTIFIC_CALIB_COEFFICIENT');
netcdf.putVar(fid,varidout,scientific_calib_coefficient);
clear scientific_calib_coefficient varidin
%
% SCIENTIFIC_CALIB_COMMENT (63)
%
varidin=netcdf.inqVarID(fidinfile,'SCIENTIFIC_CALIB_COMMENT');
scientific_calib_comment=netcdf.getVar(fidinfile,varidin); 
text256=repmat(' ',256,1);
text256(1:length(temp_scientific_com_dmqc))=temp_scientific_com_dmqc;
scientific_calib_comment(:,it,ncalib,in)=text256;
%
text256=repmat(' ',256,1);
text256(1:length(psal_scientific_com_dmqc))=psal_scientific_com_dmqc;
scientific_calib_comment(:,is,ncalib,in)=text256;
%
text256=repmat(' ',256,1);
text256(1:length(pres_scientific_com_dmqc))=pres_scientific_com_dmqc;
scientific_calib_comment(:,ip,ncalib,in)=text256;
varidout=netcdf.inqVarID(fid,'SCIENTIFIC_CALIB_COMMENT');
netcdf.putVar(fid,varidout,scientific_calib_comment);
%
% SCIENTIFIC_CALIB_DATE (64)
%
% dtext=datestr(now,30);
varidin=netcdf.inqVarID(fidinfile,'SCIENTIFIC_CALIB_DATE');
scientific_calib_date=netcdf.getVar(fidinfile,varidin); 
datetext=[dtext(1:8) dtext(10:15)];
scientific_calib_date(:,it,ncalib,in)=datetext;
scientific_calib_date(:,is,ncalib,in)=datetext;
scientific_calib_date(:,ip,ncalib,in)=datetext;
varidout=netcdf.inqVarID(fid,'SCIENTIFIC_CALIB_DATE');
netcdf.putVar(fid,varidout,scientific_calib_date);
%
% now prepare text for the history fields
%
% text fields to be used according to reference tables 
%
text_histinstitution='GE  '; % see reference table 4 in user manual
text_histaction='IP  ';      % see reference table 7 in user manual
text_histstep='ARSQ';        % see reference table 12 in user manual
text_histsoftware='OWC ';    % locally defined
text_histsoftwarerelease='2.0 '; % locally defined
text_histreference='ARGO CTD ref. database: CTD_for_DMQC_2019V02 + ARGO climatology ';  % please update whenever a new version of ref data is available
text_histdate=datetext;
text_histparameter='PSAL            '; % see note on page 103 in user manual, should actually be fill value
histstartpres=presanf;                 % should actually be fill value
histendpres=presend;                   % should actually be fill value
text_histqctest='                ';    % is fill value
histpreviousvalue=99999;               % is fill value
%
% HISTORY_INSTITUTION (48)
%
% es scheint das bei der ersten belegung eines unlimited feldes start und
% count Felder braucht ansonsten bekommt man eine fehlermeldung
% The number of input elements does not match the variable size.
%
varidin=netcdf.inqVarID(fidinfile,'HISTORY_INSTITUTION');
history_institution=netcdf.getVar(fidinfile,varidin);
nhistnew=nhistold+1;
history_institution_new=history_institution;
history_institution_new(:,in,nhistnew)=text_histinstitution;
varidout=netcdf.inqVarID(fid,'HISTORY_INSTITUTION');
netcdf.putVar(fid,varidout,[0,0,0],[4,nprof,nhistnew],history_institution_new);
%
% HISTORY_STEP (49)
%
varidin=netcdf.inqVarID(fidinfile,'HISTORY_STEP');
history_step=netcdf.getVar(fidinfile,varidin);
history_step_new=history_step;
history_step_new(:,in,nhistnew)=text_histstep;
varidout=netcdf.inqVarID(fid,'HISTORY_STEP');
netcdf.putVar(fid,varidout,history_step_new);
%
% HISTORY_SOFTWARE (50)
%
varidin=netcdf.inqVarID(fidinfile,'HISTORY_SOFTWARE');
history_software=netcdf.getVar(fidinfile,varidin);
history_software_new=history_software;
history_software_new(:,in,nhistnew)=text_histsoftware;
varidout=netcdf.inqVarID(fid,'HISTORY_SOFTWARE');
netcdf.putVar(fid,varidout,history_software_new);
%
% HISTORY_SOFTWARE_RELEASE (51)
%
varidin=netcdf.inqVarID(fidinfile,'HISTORY_SOFTWARE_RELEASE');
history_software_release=netcdf.getVar(fidinfile,varidin);
history_software_release_new=history_software_release;
history_software_release_new(:,in,nhistnew)=text_histsoftwarerelease;
varidout=netcdf.inqVarID(fid,'HISTORY_SOFTWARE_RELEASE');
netcdf.putVar(fid,varidout,history_software_release_new);
%
% HISTORY_REFERENCE (52)
%
varidin=netcdf.inqVarID(fidinfile,'HISTORY_REFERENCE');
history_reference=netcdf.getVar(fidinfile,varidin);
history_reference_new=history_reference;
history_reference_new(:,in,nhistnew)=text_histreference;
varidout=netcdf.inqVarID(fid,'HISTORY_REFERENCE');
netcdf.putVar(fid,varidout,history_reference_new);
%
% HISTORY_DATE (53)
%
varidin=netcdf.inqVarID(fidinfile,'HISTORY_DATE');
history_date=netcdf.getVar(fidinfile,varidin);
history_date_new=history_date;
history_date_new(:,in,nhistnew)=text_histdate;
varidout=netcdf.inqVarID(fid,'HISTORY_DATE');
netcdf.putVar(fid,varidout,history_date_new);
%
% HISTORY_ACTION (54)
%
varidin=netcdf.inqVarID(fidinfile,'HISTORY_ACTION');
history_action=netcdf.getVar(fidinfile,varidin);
history_action_new=history_action;
history_action_new(:,in,nhistnew)=text_histaction;
varidout=netcdf.inqVarID(fid,'HISTORY_ACTION');
netcdf.putVar(fid,varidout,history_action_new);
%
% HISTORY_PARAMETER (55)
%
varidin=netcdf.inqVarID(fidinfile,'HISTORY_PARAMETER');
history_parameter=netcdf.getVar(fidinfile,varidin);
history_parameter_new=history_parameter;
history_parameter_new(:,in,nhistnew)=text_histparameter;
varidout=netcdf.inqVarID(fid,'HISTORY_PARAMETER');
netcdf.putVar(fid,varidout,history_parameter_new);
%
% HISTORY_START_PRES (56)
%
varidin=netcdf.inqVarID(fidinfile,'HISTORY_START_PRES');
history_start_pres=netcdf.getVar(fidinfile,varidin);
history_start_pres_new=history_start_pres;
history_start_pres_new(in,nhistnew)=presanf;
if nprof>1
    history_start_pres_new(2,nhistnew)=99999;
end
varidout=netcdf.inqVarID(fid,'HISTORY_START_PRES');
netcdf.putVar(fid,varidout,history_start_pres_new);
%
% HISTORY_STOP_PRES (57)
%
varidin=netcdf.inqVarID(fidinfile,'HISTORY_STOP_PRES');
history_stop_pres=netcdf.getVar(fidinfile,varidin);
history_stop_pres_new=history_stop_pres;
history_stop_pres_new(in,nhistnew)=presend;
if nprof>1
        history_stop_pres_new(2,nhistnew)=99999;
end
varidout=netcdf.inqVarID(fid,'HISTORY_STOP_PRES');
netcdf.putVar(fid,varidout,history_stop_pres_new);
%
% HISTORY_PREVIOUS_VALUS (58)
%
varidin=netcdf.inqVarID(fidinfile,'HISTORY_PREVIOUS_VALUE');
history_previous_value=netcdf.getVar(fidinfile,varidin);
history_previous_value_new=history_previous_value;
history_previous_value_new(in,nhistnew)=histpreviousvalue;
if nprof>1
    history_previous_value_new(2,nhistnew)=99999;
end
varidout=netcdf.inqVarID(fid,'HISTORY_PREVIOUS_VALUE');
netcdf.putVar(fid,varidout,history_previous_value_new);
%
% HISTORY_QCTEST (59)
%
varidin=netcdf.inqVarID(fidinfile,'HISTORY_QCTEST');
history_qctest=netcdf.getVar(fidinfile,varidin);
history_qctest_new=history_qctest;
history_qctest_new(:,in,nhistnew)=text_histqctest;
varidout=netcdf.inqVarID(fid,'HISTORY_QCTEST');
netcdf.putVar(fid,varidout,history_qctest_new);
%
% % finished, now close files
% %
%
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
netcdf.close(fidinfile)
clear nlevel nprof ncalib nparam nhistold nhistnew
clear dtext datetext dateold
end
