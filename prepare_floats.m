% PREPARE_FLOATS loads float nc-files from your download folder, checks
% and prepares float data for ow_calibration, as well as adding flags
% and other DMQC parameters.
% by J. Even Ø. Nilsen, Ingrid Angel, Birgit Klein, and Kjell Arne Mork.
% DMQC-fun v0.9.3, jan.even.oeie.nilsen@hi.no.

clear all; close all
init_dmqc; % Paths and filenames.
grey=[.5 .5 .5]; % Colour setting for maps.

for I=1:length(download_dir)	% Loop floats

  close all
  clear hP hPqco hPqcn LO LA scientific_calib
  
   % Go to that float's working dir:
  cd([my_working_dir,'DMQC',filesep,float_names{I}]);
  
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
  ncread(trajfiles{I},'REPRESENTATIVE_PARK_PRESSURE'); snippet([int2str(round(max(ans(end))/100)*100),' m'],'park-depth');
  %ncread(proffiles{I},'PRES'); snippet([int2str(round(max(ans(:,end))/100)*100),' m'],'profile-depth');
  ncread(proffiles{I},'PRES'); 
  profdepth=round(nanmax(ans,[],'all')/100)*100; snippet([int2str(profdepth),' m'],'profile-depth');
  ncread(proffiles{I},'JULD'); t=ans; snippet([int2str(round(nanmean(diff(t,1,1)))),' days'],'cycle-time');
  ncread(metafiles{I},'DEPLOYMENT_PLATFORM'); snippet(ans(:,1)','ship','_deblank');
  ncread(metafiles{I},'PI_NAME'); snippet(ans(:,1)','PI','deblank');
  ncread(metafiles{I},'END_MISSION_STATUS'); 
  switch ans'
   case ' ', 'Active';
   case 'T', 'No more transmission received';
   case 'R', 'Retrieved';
  end
  snippet(ans,'float-status','_deblank');
  snippet([num2str(diff(dyear(t([1 end]))),'%4.2f'),' yrs'],'age');
  inCYCLE_NUMBER=ncread(infiles{I},'CYCLE_NUMBER')'; snippet(max(inCYCLE_NUMBER),'last-cycle');
  snippet(inCYCLE_NUMBER(find(diff(inCYCLE_NUMBER,1,2)>1))+1,'missing-cycles');
  % Grey list:
  [g1,g2,g3,g4,g5,g6,g7]=textread([download_parent_dir,'coriolis_greylist.csv'],'%s%s%s%s%s%s%s','delimiter',',');
  ii=find(contains(g1,float_names{I}));
  glist=''; for i=1:length(ii), glist=[glist,g2{ii(i)},',',g5{ii(i)},',',g6{ii(i)},';  ']; end
  snippet(glist,'grey-list','_');
  
  
  % Get the set minimum depth in OWC:
  snippet(MAP_P_EXCLUDE);
  
  % Transfer the operator info to the report:
  snippet(dmqc_operator(1).name,'dmqc_operator_name');
  snippet(dmqc_operator(1).orcid,'dmqc_operator_orcid');
  snippet(dmqc_operator(1).institution,'dmqc_operator_institution');
  snippet(dmqc_operator(1).address,'dmqc_operator_address');
  
  % Make text definitions for LaTeX report: 
  fid=fopen('wmonumdef.tex','w'); fprintf(fid,'%s%s%s','\newcommand{\WMOnum}{',float_names{I},'}'); fclose(fid);
  fid=fopen('floatsourcedef.tex','w'); fprintf(fid,'%s%s%s','\newcommand{\floatsource}{',float_dir,'}'); fclose(fid);
  fid=fopen('floatplotsdef.tex','w'); fprintf(fid,'%s%s%s','\newcommand{\floatplots}{',plot_dir,'}'); fclose(fid);
  fid=fopen('downloaddirdef.tex','w'); fprintf(fid,'%s%s%s','\newcommand{\downloaddir}{',download_dir{I},'}'); fclose(fid);



  % ------ INGEST FLOAT FILES INTO THE OWC SYSTEM: ----------------------------------

  % close all
  
  % % Go to that float's working dir so that all the info is put there.
  % cd([my_working_dir,'DMQC',filesep,float_names{I}]);
    
  comments='';
  nowchar=replace(datestr(now,30),'T','');
  
  % Find the complete height of matrix (not ideal method, but the
  % R-files have different length profiles and we need a matrix here):
 
  PARS=ncread(infiles{I},'PARAMETER');
  npars=size(PARS,2);	
  PRES=ncread(infiles{I},'PRES');
  m=size(PRES,1)+1; % Add one to make sure every profile ends with a NaN.
  % Read from the profile-files:
  Rfiles=edir(rootdirin{I});  % Will contain both D and R files, just
                              % like on the Coriolis server.
  Downcastfile=Rfiles(contains(Rfiles,'D.nc')); % The downcast filename.	
  Rfiles=Rfiles(~contains(Rfiles,'D.nc'));	% [] Ignore the downcast (?)
  % Each column corresponds to a file in the Rfiles list, not
  % (necessarily) to cycle numbers. Hence the output in write_D will be
  % put in the correct files regardless of their names (_<cyclenumber>)
  % since the Rfiles list is ported to WRITE_D.  [This has been
  % checked.] Take care to use CYCLE_NUMBER carefully throughout when
  % needed in plots and text!
  n=length(Rfiles);
  [PRES,PSAL,TEMP]=deal(nan(m,n)); % Originals, not to be changed at all!
  [PRESqco,PSALqco,TEMPqco]=deal(repmat(' ',m,n));
  %%%[SCIENTIFIC_CALIB_EQUATION,SCIENTIFIC_CALIB_COEFFICIENT,SCIENTIFIC_CALIB_COMMENT,SCIENTIFIC_CALIB_DATE]=deal(repmat(' ',256,npars,n));
  [CYCLE_NUMBER,LONG,LAT,JULD]=deal(nan(1,n));
  [POSqco,JULDqco,DIRECTION]=deal(repmat('',1,n));
  [PARAMETER]=deal(repmat(' ',16,npars,1,n));
  for i=1:n	% loop Rfiles
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
    % ncread(Rfiles{i},'PROFILE_PSAL_QC');	PROFILE_SAL_QC(1,i)=ans(1,1);
    ncread(Rfiles{i},'CYCLE_NUMBER')';	CYCLE_NUMBER(1,i)=ans(1,1);
    ncread(Rfiles{i},'DIRECTION')';	DIRECTION(1,i)=ans(1,1);
    ncread(Rfiles{i},'POSITION_QC')';	POSqco(1,i)=ans(1,1);
    ncread(Rfiles{i},'JULD_QC')';	JULDqco(1,i)=ans(1,1);
    ncread(Rfiles{i},'PARAMETER');
    try
      PARAMETER(:,:,i)=ans(:,:,1,1);
    catch
      PARAMETER(:,:,i)=ans(:,1:npars,1,1);
      IC=ans(:,npars+1:end,1,1);
      N=size(IC,2);
      % if N>0
	ic={'NB_SAMPLE_CTD','MTIME','TEMP_CNDC','PRES_MED','TEMP_MED','PSAL_MED','TEMP_STD','PSAL_STD'};
	for j=1:N
	  if contains(IC(:,j)',ic)
	    [scientific_calib.equation.(strip(IC(:,j)'))(i,:),...
	     scientific_calib.coefficient.(strip(IC(:,j)'))(i,:),...
	     scientific_calib.comment.(strip(IC(:,j)'))(i,:)] = deal('not applicable');
	    scientific_calib.date.(strip(IC(:,j)'))(i,:) = nowchar; 
	    % Nothing more is supposed to be done to these intermediate parameters.
	  else
	    warning([Rfiles{i},' has more parameters than ',infiles{I},...
		     ' and the normal intermediate parameters',...
		     '. Only reading ',snippet(deblank(ans(:,1:npars,1,1)'),' '),...
		     ', not ',snippet(deblank(ans(:,npars+1:end,1,1)'),' '),' !']);
	  end
	end
      % else
      % 	%%a=char(strcat(ans(:,1:npars,1,1)',{', '}))'; a=a(:)';
      % 	%%b=char(strcat(ans(:,npars+1:end,1,1)',{', '}))'; b=b(:)';
      % 	warning([Rfiles{i},' has more parameters than ',infiles{I},...
      % 		 '. Only reading ',snippet(deblank(ans(:,1:npars,1,1)'),' '),...
      % 		 ', not ',snippet(deblank(ans(:,npars+1:end,1,1)'),' '),' !']);
      % 	%PARAMETER(:,:,i)=ans(:,1:npars,1,1);
      % end
    end
  end
  % Some checks:
  %if size(PARAMETER)~=size(squeeze(PARS)), error('PARAMETER mismatch!'); end
  if any(unique(inCYCLE_NUMBER)~=CYCLE_NUMBER), error('CYCLE_NUMBER mismatch!'); end
  
  % Original number of measurements in each profile (for PROFILE_<PARAMETER>_QC in write_D):
  PROFILE_PRES_N=sum(~isnan(PRES)); PROFILE_TEMP_N=sum(~isnan(TEMP)); PROFILE_PSAL_N=sum(~isnan(PSAL));
  % Just a test on how percentage is calculated:
  % PROFILE_PSAL_QC=repmat(char(32),1,size(SALqco,2));
  % sum(SALqco=='1' | SALqco=='2' | SALqco=='5' | SALqco=='8')*100./PROFILE_PSAL_N*100;
  % PROFILE_PSAL_QC(ans==0)='F', PROFILE_PSAL_QC(ans>0)='E', PROFILE_PSAL_QC(ans>=25)='D'
  % PROFILE_PSAL_QC(ans>=50)='C', PROFILE_PSAL_QC(ans>=75)='B', PROFILE_PSAL_QC(ans==100)='A'
  
  REFERENCE_DATE_TIME=ncread(Rfiles{1},'REFERENCE_DATE_TIME'); % For datenum
  % I am only reading first column of each file! Other profiles are
  % not part of DMQC.

  % figure(2); plotyy(PRES(:,1),SAL(:,1),PRES(:,1),TEMP(:,1));
  % whos LONG LAT JULD PRES SAL TEMP CYCLE_NUMBER DIRECTION REFERENCE* HISTORY*; 
  % CYCLE_NUMBER
  % DIRECTION
  % ~isnan(PRES);all(pres(ans)==PRES(ans)); % Test for equality
  % return

  % This part is obsolete, since descending profile(s) are removed by
  % filename above, but check in case '*D.nc' isn not sufficient:
  if any(DIRECTION=='D'), error('There are descending profiles!'); end

  % New (numeric) QC flags to be added to here below:
  POSqcn  = zeros(size(POSqco));
  JULDqcn = zeros(size(JULDqco));
  PRESqcn = zeros(size(PRESqco));   
  PSALqcn  = zeros(size(PSALqco));   
  TEMPqcn = zeros(size(TEMPqco));   
  
  % More info for the report (including the whole appendix):
  ~(isnan(LONG)|isnan(LAT));wmosquares=unique(findwmo(LONG(ans),LAT(ans))); snippet(wmosquares); 
  fid=fopen('appendix.tex','w');
  for i=1:length(wmosquares)
    fprintf(fid,'%s\n','\begin{figure}[ph]');
    fprintf(fid,'%s\n','\centering');
    fprintf(fid,'%s%s%s%s%u%s\n','\includegraphics[width=\textwidth]{',tardir{1},filesep,'ctd_',wmosquares(i),'.eps}');%\climatology/argo_profiles/ctd_1600.eps}');
    fprintf(fid,'%s%s%s%s%u%s\n','\includegraphics[width=\textwidth]{',tardir{3},filesep,'argo_',wmosquares(i),'.eps}');%\climatology/argo_profiles/argo_1600.eps}');
    fprintf(fid,'%s%u%s\n','\caption{Overview of reference data in WMO square ',wmosquares(i),' which is traversed by Float~\WMOnum . Upper set of graphs are for CTD reference data, and lower set is for historical ARGO data. Colouring of positions in map pairs illustrate temporal coverage and latitude, respectively. The following TS and profile plots use the latter colormap.}');
    fprintf(fid,'%s%u%s\n','\label{reference',wmosquares(i),'}');
    fprintf(fid,'%s\n','\end{figure}');
  end
  fclose(fid);

  % Initialise the new calibration information objects (for all profiles, but for now just with no size):
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
  %
  % W21 3.2. Delayed-mode procedures for JULD, LATITUDE, LONGITUDE 
  %
  % Delayed-mode operators should check that JULD in the profiles
  % are in chronological order. Erroneous or missing JULD values should
  % be replaced with another telemetered value if available, or replaced
  % with interpolated values and marked with JULD_QC = ‘8’.
  find(diff(JULD)<0)+1;
  if any(ans) 
    JULD(ans)=NaN; JULDqcn(ans)=8; JULD=efill(JULD); %flags.new.JULD=length(ans);
    comments=[comments,'Non-chronological JULD; '];
    error('Some JULD not chronological! Check!');
  end
  REFERENCE_DATE_TIME'; time=datenum([ans(1:8),'T',ans(9:14)],'yyyymmddTHHMMSS')+JULD; DATES=dyear(time);
  % JULDqcn is set now.
  
  % Sort cycle number (not sure this should be done, so just a check and error for now):
  % Double cycle numbers messes up the rest. Pick out upcasts.
  [ans,IA]=sort(CYCLE_NUMBER); % C = A(IA);
  if ~all(ans==CYCLE_NUMBER)
    error('Cycle numbers not in succession!');
    LONG=LONG(IA); LAT=LAT(IA); DATES=DATES(IA); CYCLE_NUMBER=CYCLE_NUMBER(IA); PRES=PRES(:,IA); PSAL=PSAL(:,IA); TEMP=TEMP(:,IA);
    %CYCLE_NUMBER=ans;
  end
  snippet([int2str(CYCLE_NUMBER(1)),'-',int2str(CYCLE_NUMBER(end))],'cycle-numbers');

  % [] How to know about wrong dates not being just shifted? Relation
  % to cycle number?
  
  % Profile positions in LONGITUDE, LATITUDE should be checked for
  % outliers.  Erroneous or missing LONGITUDE, LATITTUDE values should
  % be replaced with another telemetered value if available, or replaced
  % with interpolated values and marked with POSITION_QC = ‘8’.
  % Done by visual of raw map. 
  
  % Default adjusted values (after sorting):
  PRES_ADJUSTED=PRES; TEMP_ADJUSTED=TEMP; PSAL_ADJUSTED=PSAL;

  % Plot the input data raw:
  [m,n]=size(PRES);
  PROFILE_NO=1:n;
  %DENS = sw_dens(PSAL,TEMP,PRES);	% For the plotting of warning plots.
  DENS = sw_pden(PSAL_ADJUSTED,TEMP_ADJUSTED,PRES_ADJUSTED,0);	% Potential density instead
  zerow=zeros(1,n);			% Row of zeros
  plot_profiles; 
  print(gcf,'-depsc',[outfiles{I}(1:end-4),'_raw.eps']);
  % 
  % comments=[comments,'POSITION outliers; '];
  % POSqcn()=8; LONG()=interp; LAT()=interp; % Set QC flags if outliers are seen
  % Set POSITION_QC when this happens.

  

  % W21 3.3. Delayed-mode procedures for pressure
  %
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
  else
    % Ikke APEX-float, men plot tid mot overflatetrykk, sensortrykk.
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
		      'Any ''4'', ''8'', or ''9'' are shown regardless of depth, not just at every 10th.}']);
  fprintf(fid,'%s\n','\label{fig:pres}');
  fprintf(fid,'%s\n','\end{figure}');
  %end
  fclose(fid);

  % Check for negative pressures, esp. after the correction:
  if any(PRES_ADJUSTED<0,'all'), error('Negative pressures detected!'); end
  
  
  %%%%%%% PHASE 2: Visual verification of Real-time mode QC (RTQC) flags %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  J=find(any(PRESqco=='4'|PSALqco=='4'|TEMPqco=='4'));
  N=length(J);
  fid=fopen('RTQCcheck.tex','w'); % For the report if any RTQC-flags
  if N>0		% If any RTQC flags at all
    find(J>checked{I});
    if any(ans)		% Manual control of _new_profiles_ with RTQC flags.
      for j=ans:N	% Loop only new columns with RTQC flags
	iP=find(PRESqco(:,J(j))=='4');
	iT=find(TEMPqco(:,J(j))=='4');
	iS=find(PSALqco(:,J(j))=='4');
	jj=J(j)+[-3:3]; jj=jj(0<jj&jj<=n); % keep indices inside 
	~isnan(LONG(jj))&~isnan(LAT(jj)); jj=jj(ans); % Avoid NaNs (POSqco is not applied yet)
	figure; set(gcf,'position',get(0,'screensize')); clf;
	[LO{j},LA{j}]=halo(LONG(jj),LAT(jj),.2); 
	[pres,temp,sal,long,lat]=inpolygon_referencedata(LO{j},LA{j},tardir);
	aT=subplot(1,2,1); plot(temp,pres,'color',grey); grid;
	hT=line(TEMP(:,jj),PRES(:,jj)); heT=line(TEMP(:,J(j)),PRES(:,J(j))); 
	if any(iT), hEp=line(TEMP(iT,J(j)),PRES(iT,J(j))); set(hEp,'color','m','marker','s','markersize',10,'linestyle','none'); end
	if any(iP), hEp=line(TEMP(iP,J(j)),PRES(iP,J(j))); set(hEp,'color','g','marker','s','markersize',10,'linestyle','none'); end
	axis ij; xlabel temp; ylabel pres; title(['CYCLE_NUMBER ',int2str(CYCLE_NUMBER(J(j)))],'interpreter','none');
	aS=subplot(1,2,2); plot(sal,pres,'color',grey); grid;
	hS=line(PSAL(:,jj),PRES(:,jj)); 
	heS=line(PSAL(:,J(j)),PRES(:,J(j))); 
	if any(iS), hEp=line(PSAL(iS,J(j)),PRES(iS,J(j))); set(hEp,'color','m','marker','s','markersize',10,'linestyle','none'); end
	if any(iP), hEp=line(TEMP(iP,J(j)),PRES(iP,J(j))); set(hEp,'color','g','marker','s','markersize',10,'linestyle','none'); end
	axis ij; xlabel sal; ylabel pres; 
	set([hT;hS],'color','b'); set([heT;heS],'color','r','linestyle','-');
	get(aT,'position'); ax=axes('position',[ans(1)+ans(3)-ans(3)/3-.01 ans(2)+.01 ans(3)/3 ans(4)/3]); hold on
	set(ax,'visible','off'); 
	plot(LO{j},LA{j},'color',grey); 
	plot(LONG(jj),LAT(jj),'marker','.','color','b'); 
	plot(LONG(J(j)),LAT(J(j)),'marker','s','color','r'); 
	sig(['Press window control buttons to zoom etc. Turn off when finished. Press any key to approve.']);
	set(gcf,'WindowKeyPressFcn','set(gcf,''name'',''checked'')');
	waitfor(gcf,'name');
	% After each first run after download and manual checking is
        % maybe done, update numbers for checked to last cycle in
        % INIT_DMQC, and run PREPARE_FLOATS (this script) again.
	% [√] Have to think about some smart way to set new flags
        % (also overwriting RTQC flags) in PRESqcn, TEMPqcn, and
        % PSALqcn while visually inspecting each. 
	% [√] For reversal qcn is set to 1.
	% Then, before writing, replace PREqco only where PRESqcr has
        % values. Then change PRESqco only where PRESqcn has flags.
	% Are there any cahnges to the data? Have RTQC changed the
        % data?	
	% √ This is solved in WRITE_D.
      end
      disp(['New profiles checked. Update ''checked''-variable in INIT_DMQC from ',...
	    int2str(checked{I}),' to ',int2str(n),' for Float ',float_names{I},...
	    ', and re-run prepare_floats!']);
    else		% no new profiles with RTQC flags => The multipanel plots:
      for j=1:N		% Loop all columns with RTQC flags
	iP=find(PRESqco(:,J(j))=='4');
	iT=find(TEMPqco(:,J(j))=='4');
	iS=find(PSALqco(:,J(j))=='4');
	jj=J(j)+[-3:3]; jj=jj(0<jj&jj<=n); % keep indices inside 
	~isnan(LONG(jj))&~isnan(LAT(jj)); jj=jj(ans); % Avoid NaNs (POSqco is not applied yet)
	if ismember(j,1:5:N) % Make new figure window
	  ju=j-1;		% j used in previous figure(s)
	  M=min(5,N-ju);	% rows in this figure
	  figure; shape=[1 1 1000 200*M];
	  set(gcf,'units','points','innerposition',shape,...
		  'paperunits','points','paperposition',shape,'PaperSize',shape(3:4),...
		  'PaperPositionMode','manual','RendererMode','manual','Renderer','opengl');
	end
	[LO{j},LA{j}]=halo(LONG(jj),LAT(jj),.2); 
	[pres,temp,sal,long,lat]=inpolygon_referencedata(LO{j},LA{j},tardir);
	aT=subplot(M,2,2*((j-ju)-1)+1); plot(temp,pres,'color',grey); 
	hT=line(TEMP(:,jj),PRES(:,jj)); heT=line(TEMP(:,J(j)),PRES(:,J(j))); 
	if any(iT), hEp=line(TEMP(iT,J(j)),PRES(iT,J(j))); set(hEp,'color','m','marker','s','linestyle','none'); end
	if any(iP), hEp=line(TEMP(iP,J(j)),PRES(iP,J(j))); set(hEp,'color','g','marker','s','linestyle','none'); end
	axis ij; xlabel temp; ylabel pres; title(['CYCLE_NUMBER ',int2str(CYCLE_NUMBER(J(j)))],'interpreter','none');
	aS=subplot(M,2,2*((j-ju)-1)+2); plot(sal,pres,'color',grey); 
	hS=line(PSAL(:,jj),PRES(:,jj)); 
	heS=line(PSAL(:,J(j)),PRES(:,J(j))); 
	if any(iS), hEp=line(PSAL(iS,J(j)),PRES(iS,J(j))); set(hEp,'color','m','marker','s','linestyle','none'); end
	if any(iP), hEp=line(PSAL(iP,J(j)),PRES(iP,J(j))); set(hEp,'color','g','marker','s','linestyle','none'); end
	axis ij; xlabel sal; ylabel pres; 
	set([hT;hS],'color','b'); set([heT;heS],'color','r','linestyle','-');
	get(aT,'position'); ax=axes('position',[ans(1)+ans(3)-ans(3)/3-.01 ans(2)+.01 ans(3)/3 ans(4)/3]); hold on
	set(ax,'visible','off'); 
	plot(LO{j},LA{j},'color',grey); 
	plot(LONG(jj),LAT(jj),'marker','.','color','b'); 
	plot(LONG(J(j)),LAT(J(j)),'marker','s','color','r'); 
	if ismember(j,[5:5:N N]), 
	  print(gcf,'-depsc',[outfiles{I}(1:end-4),'_RTQCcheck',int2str(ju/5+1),'.eps']); 
	end
	if j==min(5,N), % If multipanel plots are made, generate part for report with first graph:
	  fprintf(fid,'%s\n','Examples are shown in Figure~\ref{fig:RTQCcheck}. ');
	  fprintf(fid,'%s\n','\begin{figure}[H]');
	  %fprintf(fid,'%s\n','\centering');
	  fprintf(fid,'%s%u%s%u%s\n','\centerline{\includegraphics[width=1.2\textwidth,natwidth=',shape(3),',natheight=',shape(4),']{\floatsource\WMOnum_RTQCcheck1}}');
	  fprintf(fid,'%s\n','\caption{Float~\WMOnum. Profiles flagged by RTQC (red) compared to float profiles before and after (blue) and nearby reference data (grey).  Temperature in left and salinity in right panels, for which flagged data are marked with magenta squares, while pressure flags are shown as green squares. The inlays in left panels show the piece of float track shown in blue and the areas of origin of reference data used, as grey enclosures. These are the first flagged profiles for this float (or all; see Table~\protect\ref{tab:rtqc}). For geography, refer to Figure~\protect\ref{fig:float-info} Panel~1.}');
	  fprintf(fid,'%s\n','\label{fig:RTQCcheck}');
	  fprintf(fid,'%s\n','\end{figure}');
	end
      end
    end
  else		% No RTQC flags at all
    fprintf(fid,'%s\n','There were no RTQC flags to inspect from Float~\WMOnum.');
    system(['rm -f ',outfiles{I}(1:end-4),'_RTQCcheck*.eps']); % Clean up any previous plot
  end
  fclose(fid);
  if checked{I}<n
    disp(['New profiles downloaded. Update checked-variable in INIT_DMQC from ',int2str(checked{I}),' to ',int2str(n),' for Float ',float_names{I},'!']);
    snippet('This float is still active and further monitoring is required.','monitoring');
  elseif DATES(end)<dyear(now)-1/52*3 % three weeks allowance 
    ['There are no new profiles from this float since ',datestr(datenum(DATES(end),1,1),1),'.'];
    snippet(ans,'monitoring');
    disp([ans,' Has Float ',float_names{I},' been decommissioned after profile ',int2str(n),'?']);
  end 
  clear LO LA
  
  % PLOT the Hov-Möller diagrams of PTEMP and PSAL with raw data (i.e., disregarding flags)
  PTMP = sw_ptmp(PSAL,TEMP,PRES,0);
  figure(1001); set(gcf,'OuterPosition',[1 385 1026 600]);
  ccc=1:max(CYCLE_NUMBER); % All possible/expected cycle numbers
  [ppp,ttt,sss]=deal(nan(size(PRES,1),length(ccc))); % Corresponding matrices (empty)
  ppp(:,CYCLE_NUMBER)=PRES; ttt(:,CYCLE_NUMBER)=PTMP; sss(:,CYCLE_NUMBER)=PSAL; % Put data in right columns
  %[~,X,Y]=profcolor(CYCLE_NUMBER,PRES,PTMP);shading flat; 
  [~,X,Y]=profcolor(ccc,ppp,ttt);shading flat;
  X=X(:,CYCLE_NUMBER); Y=Y(:,CYCLE_NUMBER); % put coordinates back into compact format (for later marks on plot)
  axis ij; xlabel('CYCLE NUMBER'); ylabel('PRES');
  cmocean thermal; cbh=colorbar; ylabel(cbh,'PTMP');
  title(['Float ',char(float_names{I}),' Potential Temperature']);
  figure(1002); set(gcf,'OuterPosition',[1 385 1026 600]);
  %profcolor(CYCLE_NUMBER,PRES,PSAL);shading flat; 
  profcolor(ccc,ppp,sss);shading flat; 
  axis ij; xlabel('CYCLE NUMBER'); ylabel('PRES')
  cmocean haline; cbh=colorbar; ylabel(cbh,'PSAL');
  title(['Float ',char(float_names{I}),' Salinity (PSS-78)']);
  % Revisit these and PRINT them in the end.
  
  % % PLOT and PRINT the profile first-pressure series:
  % figure(1003); set(gcf,'OuterPosition',[1 385 1026 600]);
  % plot(CYCLE_NUMBER,PRES(1,:),'o'); grid;
  % xlabel('CYCLE NUMBER'); ylabel('Top-of-profile pressure  [dbar]');%ylabel('Surface pressure (PRES_0)  [dbar]');
  % title(['Top of profile pressure measurements by float ',char(float_names{I})]);
  % get(gca,'ylim'); set(gca,'ylim',[min(0,ans(1)) max(5,ans(2))]);
  % print(gcf,'-dpng',[outfiles{I}(1:end-4),'_PRES0.png']);  
  % PLOT and PRINT the profile first-pressure series:
  figure(1003); set(gcf,'OuterPosition',[1 385 1026 900]);
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
  PRESqco10=='2'; hPqco{2}=line(CYC10(ans)',PRES10(ans)'); set(hPqco{2},'color','c','marker','s','linestyle','none','tag','qcflag','userdata','RTQC flag ''2''');
  PRESqco10=='3'; hPqco{3}=line(CYC10(ans)',PRES10(ans)'); set(hPqco{3},'color','m','marker','s','linestyle','none','tag','qcflag','userdata','RTQC flag ''3''');
  PRESqco=='4'; hPqco{4}=line(CYC(ans)',PRES_ADJUSTED(ans)'); set(hPqco{4},'color','r','marker','s','linestyle','none','tag','qcflag','userdata','RTQC flag ''4''');
  PRESqco=='8'; hPqco{8}=line(CYC(ans)',PRES_ADJUSTED(ans)'); set(hPqco{8},'color','b','marker','.','linestyle','none','tag','qcflag','userdata','RTQC flag ''8''');
  PRESqco=='9'; hPqco{9}=line(CYC(ans)',PRES_ADJUSTED(ans)'); set(hPqco{9},'color','k','marker','.','linestyle','none','tag','qcflag','userdata','RTQC flag ''9''');
  title(ap(1),['Top of profile pressure measurements by float ',char(float_names{I})]);
  title(ap(2),['Every 10th pressure measurement by float ',char(float_names{I})]);
  % Legend and print of lower panel done below, after DMQC.
  
	% Apply flags from the RTQC: 
	% For reversal qcn has been set to 1. Do not reverse the RTQC flags
	% yet, because they should be shown properly in the graphics and
	% stats, but make sure reversed flags are not applied to the data here
	% (added '&<PARAMETER>qcn~=1' to this code).
	% NaN out columns:
	j=find(POSqco=='4'&POSqcn~=1 | JULDqco=='4'&JULDqcn~=1 & LAT<-90 & 90<LAT);  
	LONG(j)=NaN; LAT(j)=NaN; DATES(j)=NaN; PRES_ADJUSTED(:,j)=NaN; PSAL_ADJUSTED(:,j)=NaN; TEMP_ADJUSTED(:,j)=NaN;
	% NaN out single values (according to manual Section 3.3):
	PRESqco=='4'&PRESqcn~=1; PSAL_ADJUSTED(ans | PSALqco=='4'&PSALqcn~=1)=NaN; TEMP_ADJUSTED(ans | TEMPqco=='4'&TEMPqcn~=1)=NaN; PRES_ADJUSTED(ans)=NaN;
  

	
	
  %%%%%%% PHASE 3: Selected automated tests necessary for OWC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %DENS = sw_dens(PSAL_ADJUSTED,TEMP_ADJUSTED,PRES_ADJUSTED);	% For the plotting of warning plots.
  DENS = sw_pden(PSAL_ADJUSTED,TEMP_ADJUSTED,PRES_ADJUSTED,0);	% Potential density instead
  zerow=zeros(1,n);			% Row of zeros

  % Pressure increasing test / monotonically increasing pressure test:
  logical([zerow;diff(PRES_ADJUSTED,1,1)<=0]);  % The testvalue
  ans(PRES_ADJUSTED<=MAP_P_EXCLUDE)=logical(0); % Do not be concerned with monotonicity above this depth.
  if any(ans,'all')
    jnb=find(any(ans)); pnb=ans; nbt='Non-monotonic pressure';
    warning([nbt,'in float ',float_names{I},' cycle(s): ',int2str(CYCLE_NUMBER(jnb)),' !']);	% Only warning
    plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_pressure_warning.eps'])
    % Action: If there is a region of constant pressure, all but the
    % first of the consecutive levels of constant pressure should be
    % flagged as bad data (‘4’).
    % If there is a region where pressure reverses, all of the
    % pressures in the reversed part of the profile should be flagged
    % as bad data. 
    % All pressures flagged as bad data and all of the associated
    % temperatures and salinities should be removed.
    % NO ACTION. IMPLEMENT WHEN IT HAPPENS.
    comments=[comments,lower(nbt),' (>',int2str(MAP_P_EXCLUDE),' m); '];
    clear pnb jnb nbt
  else
    system(['rm -f ',outfiles{I}(1:end-4),'_pressure_warning.eps']); % Clean up any previous plot
  end
   
  % NEW! THE DOUBLE-POINTED SPIKE TEST:
  % Some times the spike is made up of two deviating points.
  % Test value = | (V2 + V3)/2 – (V1 + V4)/2 | – | (V4 – V1) / 2 |
  % For this test the testvalues are placed on the top point of the
  % imagined double point spike, hence, any finds will have to be
  % supplemented with a true value below as well. There is no chance
  % of looping over to next profile, since we pad with zeros at the
  % bottom of the testvalue matrix. We use the same criteria as for
  % single spikes.
  testvalue = [zerow; abs((PSAL_ADJUSTED(2:end-2,:)+PSAL_ADJUSTED(3:end-1,:))/2-(PSAL_ADJUSTED(4:end,:)+PSAL_ADJUSTED(1:end-3,:))/2) - abs((PSAL_ADJUSTED(4:end,:)-PSAL_ADJUSTED(1:end-3,:))/2) ; zerow ; zerow];
  [PRES_ADJUSTED<500 & testvalue>0.9 | PRES_ADJUSTED>=500 & testvalue>0.02 | PRES_ADJUSTED>=1000 & testvalue>0.01]; % By experience clear spikes in the Nordic Seas
  ans(find(ans)+1)=1; 
  if any(ans,'all')
    jnb=find(any(ans)); snb=ans; nbt='Double-pointed salinity spike(s)';
    warning([nbt,' in float ',float_names{I},' cycle(s): ',int2str(CYCLE_NUMBER(jnb)),' ! (removed)']);
    plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_S-double_spike_warning.eps']);
    PSAL_ADJUSTED(snb)=NaN; PSALqcn(snb)=4; 
    comments=[comments,lower(nbt),'; '];
    clear snb jnb nbt
  else
    system(['rm -f ',outfiles{I}(1:end-4),'_S-double_spike_warning.eps']);
  end
  testvalue = [zerow; abs((TEMP_ADJUSTED(2:end-2,:)+TEMP_ADJUSTED(3:end-1,:))/2-(TEMP_ADJUSTED(4:end,:)+TEMP_ADJUSTED(1:end-3,:))/2) - abs((TEMP_ADJUSTED(4:end,:)-TEMP_ADJUSTED(1:end-3,:))/2) ; zerow ; zerow];
  %[PRES_ADJUSTED<500 & testvalue>6 | PRES_ADJUSTED>=500 & testvalue>2]; % The RTQC9 limits 
  [PRES_ADJUSTED<500 & testvalue>6 | PRES_ADJUSTED>=500 & testvalue>2 | PRES_ADJUSTED>=1000 & testvalue>1]; % By experience in the Nordic Seas  
  ans(find(ans)+1)=1; 
  if any(ans,'all')
    jnb=find(any(ans)); tnb=ans; nbt='Double-pointed temperature spike(s)';
    warning([nbt,' in float ',float_names{I},' cycle(s): ',int2str(CYCLE_NUMBER(jnb)),' ! (removed)']);
    plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_T-double_spike_warning.eps']);
    TEMP_ADJUSTED(tnb)=NaN; TEMPqcn(tnb)=4;
    comments=[comments,lower(nbt),'; '];
    clear tnb jnb nbt
  else
    system(['rm -f ',outfiles{I}(1:end-4),'_T-double_spike_warning.eps']);
  end

  % Spike tests (RTQC double check and stricter DMQC check for some):
  % Test value = | V2 – (V3 + V1)/2 | – | (V3 – V1) / 2 |
  % according to EuroGOOS,  where V2 is the measurement being tested
  % as a spike, and V1 and V3 are the values above and below. 
  testvalue = [zerow; abs(PSAL_ADJUSTED(2:end-1,:)-(PSAL_ADJUSTED(3:end,:)+PSAL_ADJUSTED(1:end-2,:))/2) - abs((PSAL_ADJUSTED(3:end,:)-PSAL_ADJUSTED(1:end-2,:))/2) ; zerow];
  %[PRES_ADJUSTED<500 & testvalue>0.9 |  PRES_ADJUSTED>=500 & testvalue>0.3]; % The RTQC9 limits 
  [PRES_ADJUSTED<500 & testvalue>0.9 | PRES_ADJUSTED>=500 & testvalue>0.02 | PRES_ADJUSTED>=1000 & testvalue>0.01]; % By experience clear spikes in the Nordic Seas
  if any(ans,'all')
    jnb=find(any(ans)); snb=ans; nbt='Salinity spike(s)';
    warning([nbt,' in float ',float_names{I},' cycle(s): ',int2str(CYCLE_NUMBER(jnb)),' ! (removed)']);
    plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_S-spike_warning.eps']);
    PSAL_ADJUSTED(snb)=NaN; PSALqcn(snb)=4; 
    comments=[comments,lower(nbt),'; '];
    clear snb jnb nbt
  else
    system(['rm -f ',outfiles{I}(1:end-4),'_S-spike_warning.eps']);
  end
  testvalue = [zerow; abs(TEMP_ADJUSTED(2:end-1,:)-(TEMP_ADJUSTED(3:end,:)+TEMP_ADJUSTED(1:end-2,:))/2) - abs((TEMP_ADJUSTED(3:end,:)-TEMP_ADJUSTED(1:end-2,:))/2) ; zerow];
  %[PRES_ADJUSTED<500 & testvalue>6 | PRES_ADJUSTED>=500 & testvalue>2]; % The RTQC9 limits 
  [PRES_ADJUSTED<500 & testvalue>6 | PRES_ADJUSTED>=500 & testvalue>2 | PRES_ADJUSTED>=1000 & testvalue>1]; % By experience in the Nordic Seas  
  if any(ans,'all')
    jnb=find(any(ans)); tnb=ans; nbt='Temperature spike(s)';
    warning([nbt,' in float ',float_names{I},' cycle(s): ',int2str(CYCLE_NUMBER(jnb)),' ! (removed)']);
    plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_T-spike_warning.eps']);
    TEMP_ADJUSTED(tnb)=NaN; TEMP_ADJUSTEDqcn(tnb)=4;
    comments=[comments,lower(nbt),'; '];
    clear tnb jnb nbt
  else
    system(['rm -f ',outfiles{I}(1:end-4),'_T-spike_warning.eps']);
  end
  
  % Gradient test:
  % Test value = | V2 − (V3 + V1)/2 |
  % where V2 is the measurement being tested, and V1 and V3 are the values above and below.
  testvalue = [zerow ; abs(PSAL_ADJUSTED(2:end-1,:)-(PSAL_ADJUSTED(3:end,:)+PSAL_ADJUSTED(1:end-2,:))/2) ; zerow];
  [PRES_ADJUSTED<500 & testvalue>1.5 | PRES_ADJUSTED>=500 & testvalue>0.5]; % The RTQC9 limits 
  %[PRES_ADJUSTED<500 & testvalue>1.5 | PRES_ADJUSTED>=500 & testvalue>0.5 | PRES_ADJUSTED>=1000 & testvalue>0.02]; % By experience in the Nordic Seas  
    if any(ans,'all')
    jnb=find(any(ans));  snb=ans; nbt='Salinity gradient(s)';
    warning([nbt,' in float ',float_names{I},' cycle(s): ',int2str(CYCLE_NUMBER(jnb)),' ! (removed)']);
    plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_S-gradient_warning.eps']);
    PSAL_ADJUSTED(snb)=NaN; PSALqcn(snb)=4;
    clear snb jnb nbt
    comments=[comments,lower(nbt),'; '];
  else
    system(['rm -f ',outfiles{I}(1:end-4),'_S-gradient_warning.eps']);
  end
  testvalue = [zerow ; abs(TEMP_ADJUSTED(2:end-1,:)-(TEMP_ADJUSTED(3:end,:)+TEMP_ADJUSTED(1:end-2,:))/2) ; zerow];
  [PRES_ADJUSTED<500 & testvalue>9 |  PRES_ADJUSTED>=500 & testvalue>3]; % The RTQC9 limits 
  if any(ans,'all')
    jnb=find(any(ans)); tnb=ans; nbt='Temperature gradient(s)';
    warning([nbt,' in float ',float_names{I},' cycle(s): ',int2str(CYCLE_NUMBER(jnb)),' ! (removed)']);
    plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_T-gradient_warning.eps']);
    TEMP_ADJUSTED(tnb)=NaN; TEMP_ADJUSTEDqcn(tnb)=4;
    comments=[comments,lower(nbt),'; '];
    clear tnb jnb nbt 
  else
    system(['rm -f ',outfiles{I}(1:end-4),'_T-gradient_warning.eps']);
  end

  % Density inversion test:
  pres=PRES_ADJUSTED; pres(isnan(PSAL_ADJUSTED)|isnan(TEMP_ADJUSTED))=NaN;
  PR = pres - [zerow; nandiff(pres)/2];	% Reference pressure
  dens = sw_pden(PSAL_ADJUSTED,TEMP_ADJUSTED,pres,PR);	% Potential density
  downw=nandiff(dens)<-0.03;		% Downward test
  % DENS = sw_pden(PSAL_ADJUSTED,TEMP_ADJUSTED,pres,PR);	% Potential density
  % downw=nandiff(DENS)<-0.03;		% Downward test
  pres([downw;logical(zerow)])=NaN;	% Remove found
  PR = pres - [zerow; nandiff(pres)/2];	% Reference pressure again
  dens = sw_pden(PSAL_ADJUSTED,TEMP_ADJUSTED,pres,PR);	% Potential density again
  upw=nandiff(dens)<-0.03;		% Upward test
  logical([downw;zerow]+[zerow;upw]);	% Logical for bad data
  if any(ans,'all')
    jnb=find(any(ans)); inb=ans; nbt='Density inversion(s)';
    warning([nbt,' in float ',float_names{I},' cycle(s): ',int2str(CYCLE_NUMBER(jnb)),' !']);	% Only warning
    plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_inversion_warning.eps'])
    PSAL_ADJUSTED(inb)=NaN; TEMP_ADJUSTED(inb)=NaN; DENS(inb)=NaN; PSALqcn(inb)=4; TEMPqcn(inb)=4;
    comments=[comments,lower(nbt),'; '];
    clear inb jnb nbt
  else
    system(['rm -f ',outfiles{I}(1:end-4),'_inversion_warning.eps']);
  end

  % Now also the extra tests have resulted in flags and NaN-data.

  
  %%%%%%% PHASE 4: Visual DMQC of the variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  % Plot profiles in relation to nearby reference data (and inspect):
  [yrnow,~]=datevec(now);
  long=LONG; long>180; long(ans)=long(ans)-360;
  m_proj('albers','lon',mima(long),'lat',mima(LAT));[x,y]=m_ll2xy(long,LAT);
  % Find and put cluster numbers (i.e., groups) in sequence:
  if imag(Nclu{I}),	CLU1=kmeans([CYCLE_NUMBER',x',y'],abs(Nclu{I}))'; % Cluster also wrt time.
  else			CLU1=kmeans([x',y'],abs(Nclu{I}))';
  end
  CL=unique(CLU1,'stable'); CLU=nan(size(CLU1)); 
  for j=1:length(CL), CLU(CLU1==CL(j))=j; end 
  % Prepare the two figure windows
  figure(1004); figure(1005); shape=[1 1 800 250*abs(Nclu{I})];
  set([1004],'position',shape,'paperunits','points','paperposition',shape,'PaperOrientation','portrait','PaperPositionMode','manual');
  set([1005],'position',shape+[800 0 0 0],'paperunits','points','paperposition',shape,'PaperOrientation','portrait','PaperPositionMode','manual');
  clear a h ax ja
  for j=1:abs(Nclu{I})
    jj=find(CLU==j); 
    %tit=['Cluster ',int2str(j),' (Cycles ',int2str(min(CYCLE_NUMBER(jj))),'-',int2str(max(CYCLE_NUMBER(jj))),')'];
    tit=['Cluster ',int2str(j),' (Cycles ',zipnumstr(CYCLE_NUMBER(jj)),')'];
    [LO{j},LA{j}]=halo(LONG(jj),LAT(jj),.2); 
    [pres,temp,sal,long,lat,dates]=inpolygon_referencedata(LO{j},LA{j},tardir);
    save(['refdata_cluster',int2str(j)],'pres','temp','sal','long','lat','dates');
    figure(1004); % Profile plots
    a(j,1)=subplot(abs(Nclu{I}),2,2*(j-1)+1); plot(temp,pres,'color',grey); 
    h=line(TEMP_ADJUSTED(:,jj),PRES_ADJUSTED(:,jj)); set(h,'color','b');
    axis ij; xlabel temp; ylabel pres; title(tit);
    a(j,2)=subplot(abs(Nclu{I}),2,2*(j-1)+2); plot(sal ,pres,'color',grey);
    h=line(PSAL_ADJUSTED(:,jj),PRES_ADJUSTED(:,jj)); set(h,'color','b');
    axis ij; xlabel sal; ylabel pres; title(tit);
    figure(1005); % Trend plots
    ja(j,1)=subplot(abs(Nclu{I}),2,2*(j-1)+1); 
    g_y=yrnow-25:yrnow+1; g_d=1000:400:2000; % center around newyear/winter
    bin=bin2d(repmat(dates/1e10,size(sal,1),1),pres,sal,g_y,g_d,false);
    jbin=bin2d(repmat(dyear(datenum(1950,1,JULD(jj))),size(PSAL_ADJUSTED(:,jj),1),1),PRES_ADJUSTED(:,jj),PSAL_ADJUSTED(:,jj),g_y,g_d,false);
    errorbar(repmat(bin.x,size(bin.mean,1),1)',bin.mean',bin.var'./(bin.n'-1)); %plot(bin.x,bin.mean); 
    %errorbar(repmat(bin.x,size(bin.mean,1),1)',bin.mean',bin.var'); %plot(bin.x,bin.mean); 
    h=line(jbin.x,jbin.mean); set(h,'linestyle','none','marker','*')
    title(tit); ylabel sal; grid; xlim([g_y(1)-1,g_y(end)+1]);
    if j==1,al=legend([strcat({'Reference data '},int2str(jbin.yg(1:end-1)'),'-',int2str(jbin.yg(2:end)'),{' m'});...
		       strcat({'Float data '},int2str(jbin.yg(1:end-1)'),'-',int2str(jbin.yg(2:end)'),{' m'})]); 
    end 
    if all(isnan(jbin.mean),'all')|all(isnan(bin.mean),'all') 
      text(mean(xlim),mean(ylim),'No data to compare!','HorizontalAlignment','center','BackgroundColor',[.9 .9 .9]);
      set(al,'visible','off'); 
    else
      ylim(mima(bin.mean'+bin.var'./(bin.n'-1),bin.mean'-bin.var'./(bin.n'-1),jbin.mean));
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
  ax2=copyobj(ax,1005);						% Cluster map also on trend plot
  %get(al,'position'); set(al,'position',[0.6 ans(2:4)]);	% Move legend out of the way
  %get(ax2,'position'); set(ax2,'position',[0.6 ans(2:4)]);	% Move cluster map out of the way
  get(al,'position'); set(al,'position',[0.8-ans(3) ans(2:4)]);	% Move legend out of the way
  get(ax2,'position'); set(ax2,'position',[0.8-ans(3) ans(2:4)]);	% Move cluster map out of the way
  round(mima(get(ja,'ylim'))*100)/100; %set(ja,'ylim',ans);
  set(ja,'ytick',ans(1):.01:ans(2),'yminortick','on','yminorgrid','on');
  % Print figures to files:
  print(1004,'-depsc',[outfiles{I}(1:end-4),'_refcomp.eps'])
  print(1005,'-depsc',[outfiles{I}(1:end-4),'_trendcheck.eps'])
  % Generate figure parts for report:
  fid=fopen('refcomp.tex','w');
  %fprintf(fid,'%s\n','\begin{figure}[H]');
  %fprintf(fid,'%s\n','\centering');
  fprintf(fid,'%s%u%s%u%s\n','\centerline{\includegraphics[width=1.2\textwidth,natwidth=',shape(3),',natheight=',shape(4),']{\floatsource\WMOnum_refcomp}}');
  fprintf(fid,'%s\n','\caption{Float \WMOnum\ compared to nearby reference data. Float profiles (blue lines) are divided into clusters based on positions, and compared to nearby profiles from the reference data set (grey lines). Temperature in left and salinity in right panels, and one row per cluster. The inlay in first panel shows areas of origin of reference data used as grey enclosures around the clusters of positions on the blue float track.  For geography, refer to Figure~\protect\ref{fig:float-info} Panel~1.  }');
  fprintf(fid,'%s\n','\label{fig:refcomp}');
  %fprintf(fid,'%s\n','\end{figure}');
  fclose(fid);
  fid=fopen('trendcheck.tex','w');
  %fprintf(fid,'%s\n','\begin{figure}[H]');
  %fprintf(fid,'%s\n','\centering');
  fprintf(fid,'%s%u%s%u%s\n','\centerline{\includegraphics[width=1\textwidth,natwidth=',shape(3),',natheight=',shape(4),']{\floatsource\WMOnum_trendcheck}}');
  fprintf(fid,'%s\n','\caption{Temporal evolution of data from Float~\WMOnum\ and reference data (by cluster as in Figure~\protect\ref{fig:refcomp}). Time series of annual means of reference data in depth bins (see legend) are plotted as coloured lines with error bars representing the error of the mean (usually very small due to large ensembles), and annual bin means of float data by the same method plotted with asterisks. Annual bins are centered around new year (i.e., winter). The sketch under the legend shows areas of origin of reference data used as grey enclosures around the clusters of positions on the blue float track.  For geography, refer to Figure~\protect\ref{fig:float-info} Panel~1.  }');
  fprintf(fid,'%s\n','\label{fig:trendcheck}');
  %fprintf(fid,'%s\n','\end{figure}');
  fclose(fid);
  clear LO LA
  
  % W21: Bad data points identified by visual inspection from
  % delayed-mode analysts are recorded with PRES_ADJUSTED_QC = ‘4’ and
  % PRES_QC = '4'. For these bad data points, TEMP_QC, TEMP_ADJUSTED_QC,
  % PSAL_QC, PSAL_ADJUSTED_QC should also be set to ‘4’.  Please note
  % that whenever PARAM_ADJUSTED_QC = ‘4’, both PARAM_ADJUSTED and
  % PARAM_ADJUSTED_ERROR should be set to FillValue.
  %
  % √ This is taken care of both here for OWC and in WRITE_D.
  %
  % PRESqcn TEMPqcn PSALqcn % Change QC flags if necessary, all in this case.
  

 
  % 3.4. Delayed-mode procedures for temperature
  %
  % Subjective assessment following the same figures as for pressure (above).
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


	% Now apply new flags to the data going into OWC:
	% Reversed data is given qcn=1, so only applied above when removing
	% data based on RTQC.
	% NaN out columns:
	j=find(POSqcn==4 | JULDqcn==4 & LAT<-90 & 90<LAT);  
	LONG(j)=NaN; LAT(j)=NaN; DATES(j)=NaN; PRES_ADJUSTED(:,j)=NaN; PSAL_ADJUSTED(:,j)=NaN; TEMP_ADJUSTED(:,j)=NaN;
	% NaN out single values (according to manual, Section 3.3):
	PRESqcn==4; PSAL_ADJUSTED(ans | PSALqcn==4)=NaN;  TEMP_ADJUSTED(ans | TEMPqcn==4)=NaN; PRES_ADJUSTED(ans)=NaN;
	% At this point the variables here are to be considered _ADJUSTED variables.
	
  % 3.5. Delayed-mode procedures for salinity
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
  % Perform visual inspection as for PRES and TEMP above.
  %
  % PSALqcn % Change QC flags if necessary



  %%%%%%%%%%%%%%%% CORRECT DEEP ARGO/ARVOR FLOATS PRESSURE DEPENDENT CONDUCTIVITY BIAS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % "CPcor_new, should be used to re-compute salinity as the first step in 'D' mode."
  fid=fopen('prescondbias.tex','w');
  if strcmp(platyp,'ARVOR-D') & logical(1)
    disp(['--> Adjusting pressure dependent conductivity bias in ',ctdmodel,' on ',platyp,' float.']);
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
    ii=find(~isnan(PRES)&~isnan(PSAL)&~isnan(TEMP)&~isnan(PRES_ADJUSTED)&~isnan(TEMP_ADJUSTED));	% All parameters non-NaN because 
    %Co=gsw_C_from_SP(PSAL(ii),TEMP(ii),PRES(ii));	% gsw_ can't take any NaNs! :-(
    Coo=gsw_C_from_SP(PSAL(ii),TEMP(ii),PRES(ii));	% gsw_ can't take any NaNs! :-(
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
      disp(about);
      %capt=['Float \WMOnum. Salinity profiles from the ',num2ordinal(CYCLE_NUMBER(np)),' cycle using the manufacturer calibration ($CPcor_{SBE} =$ ',num2str(CPcor_SBE*1.0e+8,'%6.2f'),'~$\cdot 10^{-8}~dbar^{-1}$; i.e., raw data; green) and using the operator''s compressibility term and cell-gain ($CPcor_{new} =$ ',num2str(CPcor_new*1.0e+8,'%6.2f'),'~$\cdot~10^{-8}~dbar^{-1}$ ; cell-gain $=$ ',num2str((M_new-1)*1e3,'%6.2f'),'‰; magenta). The ship CTD profile used is from Cruise ',ctd_cruise,', Station ',ctd_station,' (black). Dashed horisontal line shows the minimum pressure the fit is done for.'];
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
      disp('No (valid) operator CPcor_new found.');  
      % Should have found 'ctd_*','CPcor_new','leg','np','minPRESS');
      find(max(PRES_ADJUSTED)>2500); np=ans(1); % Find first profile deep enough
      leg.float=[datestr(time(np)),', ',num2str(LONG(np),'%10.3f'),'\circE, ',num2str(LAT(np),'%10.3f'),'\circN, Argo ',float_names{I},', cycle ',int2str(np)]
      leg.ctd='no CTD';
      minPRESS=nan;
      % ctd_* are not used below in the case of no operator CPcor_new.
      CPcor_new=999; % so the below IF test fails
    end
    %     -20e-08 dbar -1 <= CPcor_new <= - 5e-08 dbar -1
    % • If CPcor_new falls outside those ranges, the DM operator should consider that the salinity
    % data are not correctable and that further investigations are required to demonstrate that
    % such CPcor_new values are valid.
    if CPcor_new <= -20e-8 | -5e-8 <= CPcor_new
      CPcor_new=CPcor_rec % use the recommended
      about=['For this float the operator has not determined, or not been able to determine, any valid CPcor\_new. ',...
	    'Instead the recommended CPcor\_new for ',ctdmodel,' is used.'];
      endtext='. ';
      disp(about);
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
    %oldleg=['Float PSAL (CPcorSBE = ' num2str(CPcor_SBE*1.0e+8,'%6.2f') ' \cdot 10^{-8} dbar^{-1})'];
    leg.par.old=[' (CPcorSBE = ',num2str(CPcor_SBE*1.0e+8,'%6.2f'),' \cdot 10^{-8} dbar^{-1})'];
    %leg.par.new=[' (CPcor\_new = ' num2str(CPcor_new*1.0e+8,'%6.2f') ' \cdot 10^{-8} dbar^{-1} ; cell-gain=',num2str((M_new-1)*1e3,'%6.2f'),'‰)'];
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
    fprintf(fid,'%s%s%s%s%s%s%s\n','For ARVOR-D floats, conductivity has to be corrected due to a too small conductivity cell compressibility term (\emph{CPcor}) in the manufacturer calibration of the SBE CTDs. ',about,' Figure~\ref{fig:prescondbias} shows the uncorrected and corrected salinities from the ',num2ordinal(CYCLE_NUMBER(np)),' cycle',endtext,' If the first cycle is not used it is probably because it is too shallow (see Figure~\ref{hovmoller}b).');
    %fprintf(fid,'%s\n','This re-computation salinity should be done as the first step in D-mode. However, to ensure the correct cycle is used for this calibration, it is done after the delayed-mode procedures for coordinates (Section~\ref{sec:DM-coordinates}). ');
    fprintf(fid,'%s\n','\begin{figure}[H]');
    fprintf(fid,'%s\n','\centering');
    fprintf(fid,'%s\n','\includegraphics[width=0.7\textwidth,natwidth=560,natheight=550]{\floatsource\WMOnum_prescondbias}');
    fprintf(fid,'%s%s%s\n','\caption{',capt,'}');
    fprintf(fid,'%s\n','\label{fig:prescondbias}');
    fprintf(fid,'%s\n','\end{figure}');
    
    % 2. Assess sensor drift and offset in 'D'mode by using PSAL_ADJUSTED_Cnew
    PSAL=PSAL_ADJUSTED;
    %PSAL=SAL;			% Save PSAL for the making of operator Cnew
				% and if no adjustment is being made
    PSAL_ADJUSTED=PSAL_ADJUSTED_Cnew;
    %SAL=PSAL_ADJUSTED_Cnew;	% For OWC and PSAL_ADJUSTED if no
                                % further adjustment is done.
    
  else
    PSAL=PSAL_ADJUSTED;
    %PSAL=SAL;			% Save PSAL for the making of operator Cnew
				% and if no adjustment is being made
    				
    PSAL_ADJUSTED_Cnew=[];	% Empty for saving and write_D
  end % deep arvor
  fclose(fid);
  

  if exist('SP') % PRES has been adjusted for surface pressure
    % Add plot to pressure figure:
    axes(ap(1));
    hSP(1)=line(SP.CN,SP.TPV); set(hSP(1),'marker','*','color','r')
    hSP(2)=line(SP.CN,SP.SP); set(hSP(2),'marker','o','color','r')
    set(hSP(1),'linestyle','none'); 
    legend([hP1 hSP],'Top of profile pressure (corrected)','Surface pressure offset','Despiked surface pressure offset','location','best');
    % Add PSAL scientific calib after potential CP_corr corrections:
    'Salinity re-calculated by using PRES_ADJUSTED and recorded in PSAL_ADJUSTED. ';
		scientific_calib.comment.PSAL(SP.IA,end+[1:size(ans,2)]) = repmat(ans,length(SP.IA),1);
  end
  
% 'PRES_ADJUSTED = PRES - dP'; 
  
  
  
  %%%%%%% FINAL PHASE: Final plots, snippets, save for OWC, and alter R-files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %DENS = sw_dens(PSAL_ADJUSTED,TEMP_ADJUSTED,PRES_ADJUSTED);	% Update density for the final plot:
  DENS = sw_pden(PSAL_ADJUSTED,TEMP_ADJUSTED,PRES_ADJUSTED,0);	% Potential density instead

  % Finally plot the clean plot.
  plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'.eps']);
  
  % Create this new variable:
  PTMP = sw_ptmp(PSAL_ADJUSTED,TEMP_ADJUSTED,PRES_ADJUSTED,0);

  %%%%%% Add flags to the Hov-Möller diagrams and the pressure plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % color:'w' variable itself
  % 	'grey' pressure
  % marker:point RTQC='2' or '3' 
  % 	square RTQC='4'
  % 	diamond reversed to '1'
  % 	circle new '4'
  figure(1001); 
  caxis(mima(PTMP));  % Set colorscale to the good data range
  find(TEMPqco=='2'); ho=line(X(ans),Y(ans)); set(ho,'marker','.','markersize',2,'color','w','linestyle','none');
  find(TEMPqco=='3'); ho=line(X(ans),Y(ans)); set(ho,'marker','.','markersize',2,'color','w','linestyle','none');
  find(TEMPqco=='4'); ho=line(X(ans),Y(ans)); set(ho,'marker','s','markersize',8,'color','w','linestyle','none');
  find(TEMPqcn==1);   hr=line(X(ans),Y(ans)); set(hr,'marker','d','markersize',8,'color','w','linestyle','none');
  find(TEMPqcn==4);   hn=line(X(ans),Y(ans)); set(hn,'marker','o','markersize',8,'color','w','linestyle','none');
  % find(PRESqco=='4'); ho=line(X(ans),Y(ans)); set(ho,'marker','s','markersize',8,'color',[.5 .5 .5],'linestyle','none');
  % find(PRESqcn==1);   hr=line(X(ans),Y(ans)); set(hr,'marker','d','markersize',8,'color',[.5 .5 .5],'linestyle','none');
  % find(PRESqcn==4);   hn=line(X(ans),Y(ans)); set(hn,'marker','o','markersize',8,'color',[.5 .5 .5],'linestyle','none');
  print(gcf,'-dpng',[outfiles{I}(1:end-4),'_HMtemp.png']);
  figure(1002); 
  caxis(mima(PSAL_ADJUSTED));
  find(PSALqco=='2');  ho=line(X(ans),Y(ans)); set(ho,'marker','.','markersize',2,'color','w','linestyle','none');
  find(PSALqco=='3');  ho=line(X(ans),Y(ans)); set(ho,'marker','.','markersize',2,'color','w','linestyle','none');
  find(PSALqco=='4');  ho=line(X(ans),Y(ans)); set(ho,'marker','s','markersize',8,'color','w','linestyle','none');
  find(PSALqcn==1);    hr=line(X(ans),Y(ans)); set(hr,'marker','d','markersize',8,'color','w','linestyle','none');
  find(PSALqcn==4);    hn=line(X(ans),Y(ans)); set(hn,'marker','o','markersize',8,'color','w','linestyle','none');
  % find(PRESqco=='4'); ho=line(X(ans),Y(ans)); set(ho,'marker','s','markersize',8,'color',[.5 .5 .5],'linestyle','none');
  % find(PRESqcn==1);   hr=line(X(ans),Y(ans)); set(hr,'marker','d','markersize',8,'color',[.5 .5 .5],'linestyle','none');
  % find(PRESqcn==4);   hn=line(X(ans),Y(ans)); set(hn,'marker','o','markersize',8,'color',[.5 .5 .5],'linestyle','none');
  print(gcf,'-dpng',[outfiles{I}(1:end-4),'_HMsal.png']);
  figure(1003); axes(ap(2)); % Go to pressure plot
  % See further above for RTQC flag markings on pressure.
  find(PRESqcn==1);    hPqcn{1}=line(X(ans),Y(ans)); set(hPqcn{1},'marker','d','markersize',8,'color','b','linestyle','none','tag','qcflag','userdata','DMQC flag ''1''');
  find(PRESqcn==4);    hPqcn{4}=line(X(ans),Y(ans)); set(hPqcn{4},'marker','o','markersize',8,'color','r','linestyle','none','tag','qcflag','userdata','DMQC flag ''4''');
  %
  try % Try, because some times there are no flags at all
    findobj(ap(2),'tag','qcflag'); legend(ans,get(ans,'userdata'),'location','best'); % Legend with any added flag marks (flexible legend using tag and userdata above)
  end
  % hfl=findobj(ap(2),'tag','qcflag')
  % legt=get(hfl,'userdata');
  % [hPqco{2},hPqco{3},hPqco{4},hPqco{8},hPqco{9},hPqcn{1},hPqcn{4}]; % flag marks
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
  %
  save(outfiles{I}',...
       'LAT','LONG','DATES','PRES','SAL','TEMP','PTMP','PROFILE_NO',...
       'JULD','CYCLE_NUMBER',...
       'PROFILE_*_N','*qco','*qcn','Rfiles','scientific_calib','PSAL_ADJUSTED_Cnew',...
       'ctdmodel','PSAL'); 
  % First row is for OWC, as well as OPERATOR_CPCOR_NEW and WRITE_D;
  % Second row is
  % Third row and PSAL is for WRITE_D;
  % Fourth row is for OPERATOR_CPCOR_NEW.
						

  
  %%%%%%% Finalise comments and make the table of flags %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Go to that float's working dir and write snippets for qc-table:
  cd([my_working_dir,'DMQC',filesep,float_names{I}]);
  if ~isempty(comments)
    comments(end-1:end)='. '; comments=['New qc=''4'' flags are due to: ',comments];  
  else
    comments='No new qc=''4'' flags were found necessary.';
  end
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
  fprintf(fid,frm,'POS ','''2''',sum(POSqco =='2','all'),sum(POSqco =='2'&POSqcn ==1,'all'),sum(POSqcn ==2,'all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(POSqco =='2'|POSqco =='2'&POSqcn ==1|POSqcn ==2,1))))));
  fprintf(fid,frm,'    ','''3''',sum(POSqco =='3','all'),sum(POSqco =='3'&POSqcn ==1,'all'),sum(POSqcn ==3,'all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(POSqco =='3'|POSqco =='3'&POSqcn ==1|POSqcn ==3,1))))));
  fprintf(fid,frm,'    ','''4''',sum(POSqco =='4','all'),sum(POSqco =='4'&POSqcn ==1,'all'),sum(POSqcn ==4,'all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(POSqco =='4'|POSqco =='4'&POSqcn ==1|POSqcn ==4,1))))));
  fprintf(fid,frm,'    ','''8''',sum(POSqco =='8','all'),sum(POSqco =='8'&POSqcn ==1,'all'),sum(POSqcn ==8,'all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(POSqco =='8'|POSqco =='8'&POSqcn ==1|POSqcn ==8,1))))));
  fprintf(fid,frm,'    ','''9''',sum(POSqco =='9','all'),sum(POSqco =='9'&POSqcn ==1,'all'),sum(POSqcn ==9,'all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(POSqco =='9'|POSqco =='9'&POSqcn ==1|POSqcn ==9,1))))));
  fprintf(fid,'%s\n','\hline'); 
  fprintf(fid,frm,'JULD','''2''',sum(JULDqco=='2','all'),sum(JULDqco=='2'&JULDqcn==1,'all'),sum(JULDqcn==2,'all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(JULDqco=='2'|JULDqco=='2'&JULDqcn==1|JULDqcn==2,1))))));
  fprintf(fid,frm,'    ','''3''',sum(JULDqco=='3','all'),sum(JULDqco=='3'&JULDqcn==1,'all'),sum(JULDqcn==3,'all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(JULDqco=='3'|JULDqco=='3'&JULDqcn==1|JULDqcn==3,1))))));
  fprintf(fid,frm,'    ','''4''',sum(JULDqco=='4','all'),sum(JULDqco=='4'&JULDqcn==1,'all'),sum(JULDqcn==4,'all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(JULDqco=='4'|JULDqco=='4'&JULDqcn==1|JULDqcn==4,1))))));
  fprintf(fid,frm,'    ','''8''',sum(JULDqco=='8','all'),sum(JULDqco=='8'&JULDqcn==1,'all'),sum(JULDqcn==8,'all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(JULDqco=='8'|JULDqco=='8'&JULDqcn==1|JULDqcn==8,1))))));
  fprintf(fid,frm,'    ','''9''',sum(JULDqco=='9','all'),sum(JULDqco=='9'&JULDqcn==1,'all'),sum(JULDqcn==9,'all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(JULDqco=='9'|JULDqco=='9'&JULDqcn==1|JULDqcn==9,1))))));
  fprintf(fid,'%s\n','\hline'); 
  fprintf(fid,frm,'PRES','''2''',sum(PRESqco=='2','all'),sum(PRESqco=='2'&PRESqcn==1,'all'),sum(PRESqcn==2,'all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(PRESqco=='2'|PRESqco=='2'&PRESqcn==1|PRESqcn==2,1))))));
  fprintf(fid,frm,'    ','''3''',sum(PRESqco=='3','all'),sum(PRESqco=='3'&PRESqcn==1,'all'),sum(PRESqcn==3,'all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(PRESqco=='3'|PRESqco=='3'&PRESqcn==1|PRESqcn==3,1))))));
  fprintf(fid,frm,'    ','''4''',sum(PRESqco=='4','all'),sum(PRESqco=='4'&PRESqcn==1,'all'),sum(PRESqcn==4,'all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(PRESqco=='4'|PRESqco=='4'&PRESqcn==1|PRESqcn==4,1))))));
  %fprintf(fid,frm,'    ','''8''',sum(PRESqco=='8','all'),sum(PRESqco=='8'&PRESqcn==1,'all'),sum(PRESqcn==8,'all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(PRESqco=='8'|PRESqco=='8'&PRESqcn==1|PRESqcn==8,1))))));
  fprintf(fid,frm,'    ','''9''',sum(PRESqco=='9','all'),sum(PRESqco=='9'&PRESqcn==1,'all'),sum(PRESqcn==9,'all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(PRESqco=='9'|PRESqco=='9'&PRESqcn==1|PRESqcn==9,1))))));
  fprintf(fid,'%s\n','\hline'); 
  fprintf(fid,frm,'PSAL ','''2''',sum(PSALqco =='2','all'),sum(PSALqco =='2'&PSALqcn ==1,'all'),sum(PSALqcn ==2,'all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(PSALqco =='2'|PSALqco =='2'&PSALqcn ==1|PSALqcn ==2,1))))));
  fprintf(fid,frm,'    ','''3''',sum(PSALqco =='3','all'),sum(PSALqco =='3'&PSALqcn ==1,'all'),sum(PSALqcn ==3,'all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(PSALqco =='3'|PSALqco =='3'&PSALqcn ==1|PSALqcn ==3,1))))));
  fprintf(fid,frm,'    ','''4''',sum(PSALqco =='4','all'),sum(PSALqco =='4'&PSALqcn ==1,'all'),sum(PSALqcn ==4,'all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(PSALqco =='4'|PSALqco =='4'&PSALqcn ==1|PSALqcn ==4,1))))));
  %fprintf(fid,frm,'    ','''8''',sum(PSALqco =='8','all'),sum(PSALqco =='8'&PSALqcn ==1,'all'),sum(PSALqcn ==8,'all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(PSALqco =='8'|PSALqco =='8'&PSALqcn ==1|PSALqcn ==8,1))))));
  fprintf(fid,frm,'    ','''9''',sum(PSALqco =='9','all'),sum(PSALqco =='9'&PSALqcn ==1,'all'),sum(PSALqcn ==9,'all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(PSALqco =='9'|PSALqco =='9'&PSALqcn ==1|PSALqcn ==9,1))))));
  fprintf(fid,'%s\n','\hline'); 
  fprintf(fid,frm,'TEMP','''2''',sum(TEMPqco=='2','all'),sum(TEMPqco=='2'&TEMPqcn==1,'all'),sum(TEMPqcn==2,'all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(TEMPqco=='2'|TEMPqco=='2'&TEMPqcn==1|TEMPqcn==2,1))))));
  fprintf(fid,frm,'    ','''3''',sum(TEMPqco=='3','all'),sum(TEMPqco=='3'&TEMPqcn==1,'all'),sum(TEMPqcn==3,'all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(TEMPqco=='3'|TEMPqco=='3'&TEMPqcn==1|TEMPqcn==3,1))))));
  fprintf(fid,frm,'    ','''4''',sum(TEMPqco=='4','all'),sum(TEMPqco=='4'&TEMPqcn==1,'all'),sum(TEMPqcn==4,'all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(TEMPqco=='4'|TEMPqco=='4'&TEMPqcn==1|TEMPqcn==4,1))))));
 % fprintf(fid,frm,'    ','''8''',sum(TEMPqco=='8','all'),sum(TEMPqco=='8'&TEMPqcn==1,'all'),sum(TEMPqcn==8,'all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(TEMPqco=='8'|TEMPqco=='8'&TEMPqcn==1|TEMPqcn==8,1))))));
  fprintf(fid,frm,'    ','''9''',sum(TEMPqco=='9','all'),sum(TEMPqco=='9'&TEMPqcn==1,'all'),sum(TEMPqcn==9,'all'),snippet(zipnumstr(CYCLE_NUMBER(find(any(TEMPqco=='9'|TEMPqco=='9'&TEMPqcn==1|TEMPqcn==9,1))))));
  fprintf(fid,'%s\n','\hline'); 
  fprintf(fid,'%s\n','\end{tabular}');
  fclose(fid);
  
  
  
end	% Loop floats





% Checklist on variables for OWC:
% LAT		√ (1×n, in decimal degrees, −ve means south of the equator, e.g. 20.5S = −20.5)
% LONG		√ (1×n, in decimal degrees, from 0 to 360, e.g. 98.5W in the eastern Pacific = 261.5E)
% DATES 	√ (1×n, in decimal year, e.g. 10 Dec 2000 = 2000.939726)
%		*Note that this date format is different from that used in the reference database.
% PRES		√ (m×n, in dbar, monotonically increasing; i.e. 1st element of the column is the
%		 shallowest pressure, and subsequent values are unique and increasing)
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

