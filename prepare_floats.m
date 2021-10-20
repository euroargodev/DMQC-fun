% This script loads float nc-files from your download folder, checks and
% prepares float data for ow_calibration, as well as adding flags and
% other DMQC parameters.

clear all; close all
init_dmqc; % Paths and filenames.
grey=[.5 .5 .5]; % Colour setting for maps.

for I=1:length(download_dir)	% Loop floats
  
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
  ncread(proffiles{I},'PRES'); snippet([int2str(round(nanmax(ans,[],'all')/100)*100),' m'],'profile-depth');
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
  ncread(infiles{I},'CYCLE_NUMBER'); snippet(max(ans),'last-cycle');
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

end


% ------ INGEST FLOAT FILES INTO THE OWC SYSTEM: ----------------------------------

for I=1:length(float_names)
  close all
  
  % Go to that float's working dir so that all the info is put there.
  cd([my_working_dir,'DMQC',filesep,float_names{I}]);
    
  comments='';

  % Find the complete height of matrix (not ideal method, but the
  % R-files have different length profiles and we need a matrix here):
  PARS=ncread(infiles{I},'PARAMETER');
  npars=size(PARS,2); 
  PRES=ncread(infiles{I},'PRES');
  m=size(PRES,1);
  % Read from the profile-files:
  Rfiles=edir(rootdirin{I});  % Will contain both D and R files, just
                              % like the Coriolis server.
  Downcastfile=Rfiles(contains(Rfiles,'D.nc')); % The downcast filename.	
  Rfiles=Rfiles(~contains(Rfiles,'D.nc'));	% [] Ignore the downcast (?)
  % Each column corresponds to a file in the Rfiles list, not
  % (necessarily) to cycle numbers. Hence the output in write_D will be
  % put in the correct files regardless of their names (_<cyclenumber>)
  % since the Rfiles list is ported to WRITE_D.  [This has been
  % checked.]
  % Take care to use CYCLE_NUMBER carefully throughout when needed in
  % plots and text.
  n=length(Rfiles);
  [PRES,SAL,TEMP]=deal(nan(m,n));
  [PRESqco,SALqco,TEMPqco]=deal(repmat(' ',m,n));
  %%%[SCIENTIFIC_CALIB_EQUATION,SCIENTIFIC_CALIB_COEFFICIENT,SCIENTIFIC_CALIB_COMMENT,SCIENTIFIC_CALIB_DATE]=deal(repmat(' ',256,npars,n));
  [PARS]=deal(repmat(' ',16,npars,1,n));
  [CYCLE_NUMBER,LONG,LAT]=deal(nan(1,n));
  [JULD]=deal(nan(n,1));
  [POSqco,JULDqco]=deal(repmat('',1,n));
  for i=1:n
    ncread(Rfiles{i},'LONGITUDE')';	LONG(1,i)=ans(1,1);
    ncread(Rfiles{i},'LATITUDE')';	LAT(1,i) =ans(1,1);
    ncread(Rfiles{i},'JULD'); 		JULD(i,1)=ans(1,1);
    ncread(Rfiles{i},'PRES');		PRES(1:size(ans,1),i)=ans(:,1);
    ncread(Rfiles{i},'PSAL');		SAL(1:size(ans,1),i)=ans(:,1);
    ncread(Rfiles{i},'TEMP');		TEMP(1:size(ans,1),i)=ans(:,1);
    ncread(Rfiles{i},'PRES_QC');	PRESqco(1:size(ans,1),i)=ans(:,1);
    ncread(Rfiles{i},'TEMP_QC');	TEMPqco(1:size(ans,1),i)=ans(:,1);
    ncread(Rfiles{i},'PSAL_QC');	SALqco(1:size(ans,1),i)=ans(:,1);
    % ncread(Rfiles{i},'PROFILE_PRES_QC');	PROFILE_PRES_QC(1,i)=ans(1,1);
    % ncread(Rfiles{i},'PROFILE_TEMP_QC');	PROFILE_TEMP_QC(1,i)=ans(1,1);
    % ncread(Rfiles{i},'PROFILE_PSAL_QC');	PROFILE_SAL_QC(1,i)=ans(1,1);
    ncread(Rfiles{i},'CYCLE_NUMBER')';	CYCLE_NUMBER(1,i)=ans(1,1);
    ncread(Rfiles{i},'DIRECTION')';	DIRECTION(1,i)=ans(1,1);
    ncread(Rfiles{i},'POSITION_QC')';	POSqco(1,i)=ans(1,1);
    ncread(Rfiles{i},'JULD_QC')';	JULDqco(1,i)=ans(1,1);
    ncread(Rfiles{i},'PARAMETER');	PARAMETER(:,:,i)=ans(:,:,1,1);
  end
  if size(PARAMETER)~=size(squeeze(PARS)), error('PARAMETER mismatch!'); end
  
  % Original number of measurements in each profile (for PROFILE_<PARAMETER>_QC in write_D):
  PROFILE_PRES_N=sum(~isnan(PRES)); PROFILE_TEMP_N=sum(~isnan(TEMP)); PROFILE_PSAL_N=sum(~isnan(SAL));
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
  %
  % % Remove all downcasts (only work on ascending profiles):
  % jA=find(DIRECTION=='A');	  % DIRECTION=='D';
  % LONG=LONG(jA);		  % LONG(ans)=NaN;
  % LAT=LAT(jA);		  	  % LAT(ans)=NaN;		  
  % JULD=JULD(jA);		  % JULD(ans)=NaN;		  
  % PRES=PRES(:,jA);		  % PRES(:,ans)=NaN;		  
  % SAL=SAL(:,jA);		  % SAL(:,ans)=NaN;		  
  % TEMP=TEMP(:,jA);		  % TEMP(:,ans)=NaN;		  
  % CYCLE_NUMBER=CYCLE_NUMBER(jA);  % CYCLE_NUMBER(ans)=NaN;	  
  % % Remove the downcasts in QC flags:
  % POSqco  = POSqco(:,jA);
  % JULDqco = JULDqco(:,jA);
  % PRESqco = PRESqco(:,jA);
  % SALqco  = SALqco(:,jA);
  % TEMPqco = TEMPqco(:,jA);
  % % jA is the column index to apply when stamping the original data

  % New (numeric) QC flags to be added to here below:
  POSqcn  = zeros(size(POSqco));
  JULDqcn = zeros(size(JULDqco));
  PRESqcn = zeros(size(PRESqco));   
  SALqcn  = zeros(size(SALqco));   
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
  scientific_calib.equation.PRES=''; scientific_calib.coefficient.PRES=''; scientific_calib.comment.PRES='';  scientific_calib.date.PRES='';
  scientific_calib.equation.TEMP=''; scientific_calib.coefficient.TEMP=''; scientific_calib.comment.TEMP='';  scientific_calib.date.TEMP='';
  scientific_calib.equation.PSAL=''; scientific_calib.coefficient.PSAL=''; scientific_calib.comment.PSAL='';  scientific_calib.date.PSAL='';
  % We know that all three will be DMQC'd now, so set the date already:
  % NO, we are not doing something to all.
  %[scientific_calib.date.PRES,scientific_calib.date.TEMP,scientific_calib.date.PSAL]=deal(replace(datestr(now,30),'T',''));
  
  
  % --------- DMQC: ----------------------------------------------------------------

  
  %%%%%%% PHASE 1: Delayed-mode procedures for coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  REFERENCE_DATE_TIME'; DATES=dyear(datenum([ans(1:8),'T',ans(9:14)],'yyyymmddTHHMMSS')+JULD');
  % JULDqcn is set now.
  
  % Sort cycle number (not sure this should be done, so just a check and error for now):
  % Double cycle numbers messes up the rest. Pick out upcasts.
  [ans,IA]=sort(CYCLE_NUMBER); % C = A(IA);
  if ~all(ans==CYCLE_NUMBER)
    error('Cycle numbers not in succession!');
    LONG=LONG(IA); LAT=LAT(IA); DATES=DATES(IA); CYCLE_NUMBER=CYCLE_NUMBER(IA); PRES=PRES(:,IA); SAL=SAL(:,IA); TEMP=TEMP(:,IA);
    %CYCLE_NUMBER=ans;
  end
  snippet([int2str(CYCLE_NUMBER(1)),'-',int2str(CYCLE_NUMBER(end))],'cycle-numbers');

  % [] How to know about wrong dates not being just shifted? Relation
  % to cycle number?
  
  % Profile positions in LONGITUDE, LATITUDE should be checked for
  % outliers.  Erroneous or missing LONGITUDE, LATITTUDE values should
  % be replaced with another telemetered value if available, or replaced
  % with interpolated values and marked with POSITION_QC = ‘8’.
  % Done by visual of raw map. Plot the input data raw:
  [m,n]=size(PRES);
  PROFILE_NO=1:n;
  %DENS = sw_dens(SAL,TEMP,PRES);	% For the plotting of warning plots.
  DENS = sw_pden(SAL,TEMP,PRES,0);	% Potential density instead
  zerow=zeros(1,n);			% Row of zeros
  plot_profiles; 
  print(gcf,'-depsc',[outfiles{I}(1:end-4),'_raw.eps']);
  % 
  % comments=[comments,'POSITION outliers; '];
  % POSqcn()=8; LONG()=interp; LAT()=interp; % Set QC flags if outliers are seen
  % Set POSITION_QC when this happens.

    
  %%%%%%% PHASE 2: Visual verification of Real-time mode QC (RTQC) flags %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  J=find(any(PRESqco=='4'|SALqco=='4'|TEMPqco=='4'));
  N=length(J);
  fid=fopen('RTQCcheck.tex','w'); % For the report if any RTQC-flags
  if N>0		% If any RTQC flags at all
    find(J>checked{I});
    if any(ans)		% Manual control of _new_profiles_ with RTQC flags.
      for j=ans:N	% Loop only new columns with RTQC flags
	iP=find(PRESqco(:,J(j))=='4');
	iT=find(TEMPqco(:,J(j))=='4');
	iS=find(SALqco(:,J(j))=='4');
	jj=J(j)+[-3:3]; jj=jj(0<jj&jj<=n);
	figure; set(gcf,'position',get(0,'screensize')); clf;
	[LO{j},LA{j}]=halo(LONG(jj),LAT(jj),.2); 
	[pres,temp,sal,long,lat]=inpolygon_referencedata(LO{j},LA{j},tardir);
	aT=subplot(1,2,1); plot(temp,pres,'color',grey); grid;
	hT=line(TEMP(:,jj),PRES(:,jj)); heT=line(TEMP(:,J(j)),PRES(:,J(j))); 
	if any(iT), hEp=line(TEMP(iT,J(j)),PRES(iT,J(j))); set(hEp,'color','m','marker','s','markersize',10,'linestyle','none'); end
	if any(iP), hEp=line(TEMP(iP,J(j)),PRES(iP,J(j))); set(hEp,'color','g','marker','s','markersize',10,'linestyle','none'); end
	axis ij; xlabel temp; ylabel pres; title(['CYCLE_NUMBER ',int2str(CYCLE_NUMBER(J(j)))],'interpreter','none');
	aS=subplot(1,2,2); plot(sal,pres,'color',grey); grid;
	hS=line(SAL(:,jj),PRES(:,jj)); 
	heS=line(SAL(:,J(j)),PRES(:,J(j))); 
	if any(iS), hEp=line(SAL(iS,J(j)),PRES(iS,J(j))); set(hEp,'color','m','marker','s','markersize',10,'linestyle','none'); end
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
        % SALqcn while visually inspecting each. 
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
	iS=find(SALqco(:,J(j))=='4');
	jj=J(j)+[-3:3]; jj=jj(0<jj&jj<=n);
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
	hS=line(SAL(:,jj),PRES(:,jj)); 
	heS=line(SAL(:,J(j)),PRES(:,J(j))); 
	if any(iS), hEp=line(SAL(iS,J(j)),PRES(iS,J(j))); set(hEp,'color','m','marker','s','linestyle','none'); end
	if any(iP), hEp=line(TEMP(iP,J(j)),PRES(iP,J(j))); set(hEp,'color','g','marker','s','linestyle','none'); end
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
	  fprintf(fid,'%s\n','\caption{Float~\WMOnum. Profiles flagged by RTQC (red) compared to float profiles before and after (blue) and nearby reference data (grey).  Temperature in left and salinity in right panels. Flagged data are marked with magenta squares. If pressure is the flagged variable, the circle is put on both the temperature and salinity profile. The inlays in left panels show the piece of float track shown in blue and the areas of origin of reference data used, as grey enclosures. These are the first flagged profiles for this float (or all; see Table~\protect\ref{tab:rtqc}). For geography, refer to Figure~\protect\ref{fig:float-info} Panel~1.}');
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
  
  % PLOT the Hov-Möller diagrams of PTEMP and PSAL with raw data (i.e., disregarding flags)
  PTMP = sw_ptmp(SAL,TEMP,PRES,0);
  figure(1001); set(gcf,'OuterPosition',[1 385 1026 600]);
  %pcolor(CYCLE_NUMBER-.5,PRES,PTMP);shading flat; 
  [~,X,Y]=profcolor(CYCLE_NUMBER,PRES,PTMP);shading flat; 
  axis ij; xlabel('CYCLE NUMBER'); ylabel('PRES');
  cmocean thermal; cbh=colorbar; ylabel(cbh,'PTMP');
  title(['Float ',char(float_names{I}),' Potential Temperature']);
  figure(1002); set(gcf,'OuterPosition',[1 385 1026 600]);
  %pcolor(CYCLE_NUMBER-.5,PRES,SAL);shading flat; 
  profcolor(CYCLE_NUMBER,PRES,SAL);shading flat; 
  axis ij; xlabel('CYCLE NUMBER'); ylabel('PRES')
  cmocean haline; cbh=colorbar; ylabel(cbh,'PSAL');
  title(['Float ',char(float_names{I}),' Salinity (PSS-78)']);
  % Revisit these and PRINT them in the end.
  
  % PLOT and PRINT the profile first-pressure series:
  figure(1003); set(gcf,'OuterPosition',[1 385 1026 600]);
  plot(CYCLE_NUMBER,PRES(1,:),'o'); grid;
  xlabel('CYCLE NUMBER'); ylabel('Top-of-profile pressure  [dbar]');%ylabel('Surface pressure (PRES_0)  [dbar]');
  title(['Top of profile pressure measurements by float ',char(float_names{I})]);
  get(gca,'ylim'); set(gca,'ylim',[min(0,ans(1)) max(5,ans(2))]);
  print(gcf,'-dpng',[outfiles{I}(1:end-4),'_PRES0.png']);  
  
  
	% Apply flags from the RTQC: 
	% For reversal qcn has been set to 1. Do not reverse the RTQC flags
	% yet, because they should be shown properly in the graphics and
	% stats, but make sure reversed flags are not applied to the data here
	% (adding '&POSqcn~=1').
	% NaN out columns:
	j=find(POSqco=='4'&POSqcn~=1 | JULDqco=='4'&JULDqcn~=1 & LAT<-90 & 90<LAT);  
	LONG(j)=NaN; LAT(j)=NaN; DATES(j)=NaN; PRES(:,j)=NaN; SAL(:,j)=NaN; TEMP(:,j)=NaN;
	% NaN out single values (according to manual Section 3.3):
	PRESqco=='4'&PRESqcn~=1; SAL(ans | SALqco=='4'&SALqcn~=1)=NaN; TEMP(ans | TEMPqco=='4'&TEMPqcn~=1)=NaN; PRES(ans)=NaN;
  
  
  %%%%%%% PHASE 3: Selected automated tests necessary for OWC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %DENS = sw_dens(SAL,TEMP,PRES);	% For the plotting of warning plots.
  DENS = sw_pden(SAL,TEMP,PRES,0);	% Potential density instead
  zerow=zeros(1,n);			% Row of zeros

  % Pressure increasing test / monotonically increasing pressure test:
  logical([zerow;diff(PRES,1,1)<=0]); 
  ans(PRES<=MAP_P_EXCLUDE)=logical(0); % Do not be concerned with monotonicity above this depth.
  if any(ans,'all')
    jnb=find(any(ans)); pnb=ans; nbt='Non-monotonic pressure';
    warning([nbt,'in float ',float_names{I},' profile(s): ',int2str(PROFILE_NO(jnb)),' !']);	% Only warning
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
  testvalue = [zerow; abs((SAL(2:end-2,:)+SAL(3:end-1,:))/2-(SAL(4:end,:)+SAL(1:end-3,:))/2) - abs((SAL(4:end,:)-SAL(1:end-3,:))/2) ; zerow ; zerow];
  [PRES<500 & testvalue>0.9 | PRES>=500 & testvalue>0.02 | PRES>=1000 & testvalue>0.01]; % By experience clear spikes in the Nordic Seas
  ans(find(ans)+1)=1; 
  if any(ans,'all')
    jnb=find(any(ans)); snb=ans; nbt='Double-pointed salinity spike(s)';
    warning([nbt,' in float ',float_names{I},' profile(s): ',int2str(PROFILE_NO(jnb)),' ! (removed)']);
    plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_S-double_spike_warning.eps']);
    SAL(snb)=NaN; SALqcn(snb)=4; 
    comments=[comments,lower(nbt),'; '];
    clear snb jnb nbt
  else
    system(['rm -f ',outfiles{I}(1:end-4),'_S-double_spike_warning.eps']);
  end
  testvalue = [zerow; abs((TEMP(2:end-2,:)+TEMP(3:end-1,:))/2-(TEMP(4:end,:)+TEMP(1:end-3,:))/2) - abs((TEMP(4:end,:)-TEMP(1:end-3,:))/2) ; zerow ; zerow];
  %[PRES<500 & testvalue>6 | PRES>=500 & testvalue>2]; % The RTQC9 limits 
  [PRES<500 & testvalue>6 | PRES>=500 & testvalue>2 | PRES>=1000 & testvalue>1]; % By experience in the Nordic Seas  
  ans(find(ans)+1)=1; 
  if any(ans,'all')
    jnb=find(any(ans)); tnb=ans; nbt='Double-pointed temperature spike(s)';
    warning([nbt,' in float ',float_names{I},' profile(s): ',int2str(PROFILE_NO(jnb)),' ! (removed)']);
    plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_T-double_spike_warning.eps']);
    TEMP(tnb)=NaN; TEMPqcn(tnb)=4;
    comments=[comments,lower(nbt),'; '];
    clear tnb jnb nbt
  else
    system(['rm -f ',outfiles{I}(1:end-4),'_T-double_spike_warning.eps']);
  end

  % Spike tests (RTQC double check and stricter DMQC check for some):
  % Test value = | V2 – (V3 + V1)/2 | – | (V3 – V1) / 2 |
  % according to EuroGOOS,  where V2 is the measurement being tested
  % as a spike, and V1 and V3 are the values above and below. 
  testvalue = [zerow; abs(SAL(2:end-1,:)-(SAL(3:end,:)+SAL(1:end-2,:))/2) - abs((SAL(3:end,:)-SAL(1:end-2,:))/2) ; zerow];
  %[PRES<500 & testvalue>0.9 |  PRES>=500 & testvalue>0.3]; % The RTQC9 limits 
  [PRES<500 & testvalue>0.9 | PRES>=500 & testvalue>0.02 | PRES>=1000 & testvalue>0.01]; % By experience clear spikes in the Nordic Seas
  if any(ans,'all')
    jnb=find(any(ans)); snb=ans; nbt='Salinity spike(s)';
    warning([nbt,' in float ',float_names{I},' profile(s): ',int2str(PROFILE_NO(jnb)),' ! (removed)']);
    plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_S-spike_warning.eps']);
    SAL(snb)=NaN; SALqcn(snb)=4; 
    comments=[comments,lower(nbt),'; '];
    clear snb jnb nbt
  else
    system(['rm -f ',outfiles{I}(1:end-4),'_S-spike_warning.eps']);
  end
  testvalue = [zerow; abs(TEMP(2:end-1,:)-(TEMP(3:end,:)+TEMP(1:end-2,:))/2) - abs((TEMP(3:end,:)-TEMP(1:end-2,:))/2) ; zerow];
  %[PRES<500 & testvalue>6 | PRES>=500 & testvalue>2]; % The RTQC9 limits 
  [PRES<500 & testvalue>6 | PRES>=500 & testvalue>2 | PRES>=1000 & testvalue>1]; % By experience in the Nordic Seas  
  if any(ans,'all')
    jnb=find(any(ans)); tnb=ans; nbt='Temperature spike(s)';
    warning([nbt,' in float ',float_names{I},' profile(s): ',int2str(PROFILE_NO(jnb)),' ! (removed)']);
    plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_T-spike_warning.eps']);
    TEMP(tnb)=NaN; TEMPqcn(tnb)=4;
    comments=[comments,lower(nbt),'; '];
    clear tnb jnb nbt
  else
    system(['rm -f ',outfiles{I}(1:end-4),'_T-spike_warning.eps']);
  end
  
  % Gradient test:
  % Test value = | V2 − (V3 + V1)/2 |
  % where V2 is the measurement being tested, and V1 and V3 are the values above and below.
  testvalue = [zerow ; abs(SAL(2:end-1,:)-(SAL(3:end,:)+SAL(1:end-2,:))/2) ; zerow];
  [PRES<500 & testvalue>1.5 | PRES>=500 & testvalue>0.5]; % The RTQC9 limits 
  %[PRES<500 & testvalue>1.5 | PRES>=500 & testvalue>0.5 | PRES>=1000 & testvalue>0.02]; % By experience in the Nordic Seas  
    if any(ans,'all')
    jnb=find(any(ans));  snb=ans; nbt='Salinity gradient(s)';
    warning([nbt,' in float ',float_names{I},' profile(s): ',int2str(PROFILE_NO(jnb)),' ! (removed)']);
    plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_S-gradient_warning.eps']);
    SAL(snb)=NaN; SALqcn(snb)=4;
    clear snb jnb nbt
    comments=[comments,lower(nbt),'; '];
  else
    system(['rm -f ',outfiles{I}(1:end-4),'_S-gradient_warning.eps']);
  end
  testvalue = [zerow ; abs(TEMP(2:end-1,:)-(TEMP(3:end,:)+TEMP(1:end-2,:))/2) ; zerow];
  [PRES<500 & testvalue>9 |  PRES>=500 & testvalue>3]; % The RTQC9 limits 
  if any(ans,'all')
    jnb=find(any(ans)); tnb=ans; nbt='Temperature gradient(s)';
    warning([nbt,' in float ',float_names{I},' profile(s): ',int2str(PROFILE_NO(jnb)),' ! (removed)']);
    plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_T-gradient_warning.eps']);
    TEMP(tnb)=NaN; TEMPqcn(tnb)=4;
    comments=[comments,lower(nbt),'; '];
    clear tnb jnb nbt 
  else
    system(['rm -f ',outfiles{I}(1:end-4),'_T-gradient_warning.eps']);
  end

  % Density inversion test:
  pres=PRES; pres(isnan(SAL)|isnan(TEMP))=NaN;
  PR = pres - [zerow; nandiff(pres)/2];	% Reference pressure
  dens = sw_pden(SAL,TEMP,pres,PR);	% Potential density
  downw=nandiff(dens)<-0.03;		% Downward test
  % DENS = sw_pden(SAL,TEMP,pres,PR);	% Potential density
  % downw=nandiff(DENS)<-0.03;		% Downward test
  pres([downw;logical(zerow)])=NaN;	% Remove found
  PR = pres - [zerow; nandiff(pres)/2];	% Reference pressure again
  dens = sw_pden(SAL,TEMP,pres,PR);	% Potential density again
  upw=nandiff(dens)<-0.03;		% Upward test
  logical([downw;zerow]+[zerow;upw]);	% Logical for bad data
  if any(ans,'all')
    jnb=find(any(ans)); inb=ans; nbt='Density inversion(s)';
    warning([nbt,' in float ',float_names{I},' profile(s): ',int2str(PROFILE_NO(jnb)),' !']);	% Only warning
    plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_inversion_warning.eps'])
    SAL(inb)=NaN; TEMP(inb)=NaN; DENS(inb)=NaN; SALqcn(inb)=4; TEMPqcn(inb)=4;
    comments=[comments,lower(nbt),'; '];
    clear inb jnb nbt
  else
    system(['rm -f ',outfiles{I}(1:end-4),'_inversion_warning.eps']);
  end

  % Now also the extra tests have resulted in flags and NaN-data.



  %%%%%%% PHASE 4: Visual DMQC of the variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % W21 3.3. Delayed-mode procedures for pressure
  %
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
    h=line(TEMP(:,jj),PRES(:,jj)); set(h,'color','b');
    axis ij; xlabel temp; ylabel pres; title(tit);
    a(j,2)=subplot(abs(Nclu{I}),2,2*(j-1)+2); plot(sal ,pres,'color',grey);
    h=line(SAL(:,jj),PRES(:,jj)); set(h,'color','b');
    axis ij; xlabel sal; ylabel pres; title(tit);
    figure(1005); % Trend plots
    g_y=yrnow-25:yrnow; g_d=1000:400:2000; % center around newyear/winter
    bin=bin2d(repmat(dates/1e10,size(sal,1),1),pres,sal,g_y,g_d,false);
    jbin=bin2d(repmat(dyear(datenum(1950,1,JULD(jj)))',size(SAL(:,jj),1),1),PRES(:,jj),SAL(:,jj),g_y,g_d,false);
    ja(j,1)=subplot(abs(Nclu{I}),2,2*(j-1)+1); 
    errorbar(repmat(bin.x,size(bin.mean,1),1)',bin.mean',bin.var'./(bin.n'-1)); %plot(bin.x,bin.mean); 
    %errorbar(repmat(bin.x,size(bin.mean,1),1)',bin.mean',bin.var'); %plot(bin.x,bin.mean); 
    h=line(jbin.x,jbin.mean); set(h,'linestyle','none','marker','*')
    title(tit); ylabel sal; grid; xlim([g_y(1)-1,g_y(end)+1]);
    if j==1,al=legend([strcat({'Reference data '},int2str(jbin.yg(1:end-1)'),'-',int2str(jbin.yg(2:end)'),{' m'});...
		       strcat({'Float data '},int2str(jbin.yg(1:end-1)'),'-',int2str(jbin.yg(2:end)'),{' m'})]); end 
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
  get(al,'position'); set(al,'position',[0.6 ans(2:4)]);	% Move legend out of the way
  get(ax2,'position'); set(ax2,'position',[0.6 ans(2:4)]);	% Move cluster map out of the way
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
  
  % W21: Bad data points identified by visual inspection from
  % delayed-mode analysts are recorded with PRES_ADJUSTED_QC = ‘4’ and
  % PRES_QC = '4'. For these bad data points, TEMP_QC, TEMP_ADJUSTED_QC,
  % PSAL_QC, PSAL_ADJUSTED_QC should also be set to ‘4’.  Please note
  % that whenever PARAM_ADJUSTED_QC = ‘4’, both PARAM_ADJUSTED and
  % PARAM_ADJUSTED_ERROR should be set to FillValue.
  %
  % √ This is taken care of both here for OWC and in WRITE_D.
  %
  % PRESqcn TEMPqcn SALqcn % Change QC flags if necessary, all in this case.
  
  % Pressure adjustments for APEX, or not:
  fid=fopen('pressure-adjustment.tex','w');
  if contains(platyp,'APEX') 
    fprintf(fid,'%s%s%s\n','The pressure adjustments for ',...
	    platyp,...
	    ' floats (Wong et al., 2021) \textbf{are not implemented yet}. ');

    % 3.3.1. Delayed-mode pressure adjustment for APEX floats
    % try,   SP=ncread(techfiles{I},'PRES_SurfaceOffsetNotTruncated_dbar'); % Surface pressure
    % catch, SP=ncread(techfiles{I},'PRES_SurfaceOffsetTruncatedPlus5dbar_dbar'); end  % Surface pressure
    % error('Please implement delayed-mode pressure adjustment for APEX floats!');
  
    % 3.3.2. Truncated negative pressure drift (TNPD) in APEX floats
    % error('Please implement truncated negative pressure drift (TNPD) in APEX floats!');

    % scientific_calib.equation.PRES	= repmat('PRES_ADJUSTED = PRES – dP',n,1); % First new entry
    % scientific_calib.coefficient.PRES	= 'dP = -0.2 dbar'; % Maybe different between profiles?
    % scientific_calib.comment.PRES	= repmat('Pressures adjusted by using pressure offset at the sea surface. The quoted error is manufacturer specified accuracy in dbar.',n,1); % First new entry
    % scientific_calib.date.PRES	= repmat(replace(datestr(now,30),'T',''),n,1); % First new entry
    
  else
    % Ikke APEX-float, men plot tid mot overflatetrykk, sensortrykk.
    fprintf(fid,'%s%s%s\n','Sea surface pressure adjustments should be done for APEX floats (Wong et al., 2021). ');
    fprintf(fid,'%s%s%s\n','This is an ',...
	    platyp,...
	    ' float. Instead, Figure~\ref{surf_press} shows the surface pressure from the top of each profile. ');
    
    % Assuming that scientific_calib should be 'none' throughout:
    scientific_calib.equation.PRES	= repmat('none',n,1); % First new entry 
    scientific_calib.coefficient.PRES	= repmat('none',n,1); % First new entry
    scientific_calib.comment.PRES	= repmat('none',n,1); % First new entry
    scientific_calib.date.PRES		= repmat(replace(datestr(now,30),'T',''),n,1); % First new entry
  
  end

  fprintf(fid,'%s\n','\begin{figure}[H]');
  fprintf(fid,'%s\n','\centering');
  fprintf(fid,'%s\n','\includegraphics[width=\textwidth,natwidth=1500,natheight=750]{\floatsource\WMOnum_PRES0.png}');
  fprintf(fid,'%s\n','\caption{Float \WMOnum. Top of profile pressure series. Blue circles indicate pressure value in the real-time.}');
  %    \includegraphics[width=\textwidth]{Example_float/surf_pres_\WMOnum}
  %    \caption{Float \WMOnum. Sea surface pressure data. The red crosses indicate the raw pressure before float descent, recorded after sending data to GDAC. Blue circle indicate pressure value in the real-time. Green rotated cross shows the pressure correction applied from the previous float cycle. Top plot- data constrained between -2.4 and 2.4 dbar, middle plot- data constrained between -20 and +20 dbar, bottom plot- data with a max range of data.}
  fprintf(fid,'%s\n','\label{surf_press}');
  fprintf(fid,'%s\n','\end{figure}');
  %end
  fclose(fid);
  % [] Check if this needs to be done earlier or later!


 
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
  % When data are good and no adjustment is needed:
  scientific_calib.equation.TEMP	= repmat('none',n,1); % First new entry 
  scientific_calib.coefficient.TEMP	= repmat('none',n,1); % First new entry
  scientific_calib.comment.TEMP		= repmat('The quoted error is manufacturer specified accuracy with respect to ITS-90 at time of laboratory calibration.',n,1); % First new entry
  scientific_calib.date.TEMP		= repmat(replace(datestr(now,30),'T',''),n,1); % First new entry
  %
  % TEMPqc % Change QC flags if necessary


	% Now apply new flags to the data going into OWC:
	% Reversed data is given qcn=1, so only applied above when removing
	% data based on RTQC.
	% NaN out columns:
	j=find(POSqcn==4 | JULDqcn==4 & LAT<-90 & 90<LAT);  
	LONG(j)=NaN; LAT(j)=NaN; DATES(j)=NaN; PRES(:,j)=NaN; SAL(:,j)=NaN; TEMP(:,j)=NaN;
	% NaN out single values (according to manual, Section 3.3):
	PRESqcn==4; SAL(ans | SALqcn==4)=NaN;  TEMP(ans | TEMPqcn==4)=NaN; PRES(ans)=NaN;
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
  % SALqcn % Change QC flags if necessary



  %%%%%%%%%%%%%%%% CORRECT DEEP ARGO/ARVOR FLOATS PRESSURE DEPENDENT CONDUCTIVITY BIAS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % "CPcor_new, should be used to re-compute salinity as the first step in 'D' mode."
  fid=fopen('prescondbias.tex','w');
  if strcmp(platyp,'ARVOR-D') & logical(1)
    disp(['--> Adjusting pressure dependent conductivity bias in ',ctdmodel,' on ',platyp,' float.']);
    % The steps to re-compute salinity with CPcor_new are as follows.
    % (a). Fill PRES_ADJUSTED and TEMP_ADJUSTED [] ??????
    % Follow the same 'D' mode procedures as for the 2000-dbar floats,
    % as in Sections 3.3 and 3.4.  
    % [√] Should be PRES and TEMP as we now have here. 
    % (b). Compute original conductivity, Co.  This is done by using
    % PRES, TEMP, PSAL. For example, if using the Gibbs-SeaWater (GSW)
    % Oceanographic Matlab Toolbox, then Co = gsw_C_from_SP (PSAL, TEMP,
    % PRES).
    ii=find(~isnan(PRES)&~isnan(SAL)&~isnan(TEMP));	% All parameters non-NaN because 
    Coo=gsw_C_from_SP(SAL(ii),TEMP(ii),PRES(ii));	% gsw_ can't take any NaNs! :-(
    Co=nan(m,n); Co(ii)=Coo; clear Coo ii		% Matrix again
    % (c). Compute new conductivity, Cnew.
    % Cnew = Co*(1 + d*TEMP + CPcor_SBE*PRES) / (1 + d*TEMP_ADJUSTED +
    % CPcor_new*PRES_ADJUSTED); where
    % d = 3.25e-06,
    % CPcor_SBE = -9.57e-08 dbar-1,
    % CPcor_new = -12.5e-08 dbar-1 for SBE-61 data,
    % CPcor_new = -13.5e-08 dbar-1 for Deep SBE-41CP data, or
    % CPcor_new = as estimated by DMQC operator.
    CPcor_SBE = -9.57e-8;
    switch ctdmodel
     case 'SBE41CP', CPcor_new = -13.5e-8;
     case 'SBE61',  CPcor_new = -12.5e-8;
    end
    % Or CPcor_new = as estimated by DMQC operator.
    d = 3.25e-6;
    Cnew = Co.*(1 + d*TEMP + CPcor_SBE*PRES) ./ (1 + d*TEMP + CPcor_new*PRES); 
    % (d). Compute new salinity, PSAL_ADJUSTED_Cnew.  This is done by
    % using PRES_ADJUSTED, TEMP_ADJUSTED, and Cnew. For example, if
    % using the Gibbs-SeaWater (GSW) Oceanographic Matlab Toolbox, then
    % PSAL_ADJUSTED_Cnew = gsw_SP_from_C (Cnew, TEMP_ADJUSTED,
    % PRES_ADJUSTED).
    PSAL_ADJUSTED_Cnew=gsw_SP_from_C(Cnew,TEMP,PRES);
    % [] What about cell-gain? What is it? 
    figure(2001);clf;set(gcf,'OuterPosition',[1027 385 400 800]);	% Look at correction of first profile.  
    j=1; 
    plot(PRES(:,j),SAL(:,j),'.',PRES(:,j),PSAL_ADJUSTED_Cnew(:,j),'.'); 
    % i=find(PRES(:,j)>=500);
    % plot(PRES(i,j),SAL(i,j),'.',PRES(i,j),PSAL_ADJUSTED_Cnew(i,j),'.'); 
    % mima(SAL(i(end),j),PSAL_ADJUSTED_Cnew(i(end),j)); ylim(ans+diff(ans)*[-1 1]);
    % xlim;xlim([500 ans(2)]);
    ylim([34.895 34.925]); xlim([0 3500]);
    grid; xlabel PRES; ylabel S;  view([90 90]);
    legend(['PSAL'],...
	   ['PSAL_ADJUSTED_Cnew'],...
	   'interpreter','none','location','best'); 
    title(['Pressure dependent conductivity-bias correction (Cycle ',int2str(CYCLE_NUMBER(j)),')']);
    print(gcf,'-depsc',[outfiles{I}(1:end-4),'_prescondbias.eps']);
    % figure(2002);clf;	% Look at the stability of differences (maybe void)
    % find(all(~isnan(PSAL_ADJUSTED_Cnew),2));i=ans(end-200:end); bias=PSAL_ADJUSTED_Cnew(i,:)-SAL(i,:); boxplot(bias);
    % mima(PRES(i,:)); title(['PSAL_ADJUSTED_Cnew - original PSAL (in depth range ',int2str(ans(1)),'-',int2str(ans(2)),' m)'],'interpreter','none');
    % We really need CTD refereces in order to chack anything.
    % [] Is there something to learn from assumning stable sal or dens in the deep layers?
    
    PSAL_ADJUSTED=SAL;		% Save SAL_ADJUSTED before making the Cnew one.
    SAL=PSAL_ADJUSTED_Cnew;	% For OWC.
    
    % Calibration info:
    scientific_calib_equation.PSAL	= repmat('New conductivity = original conductivity * (1 + delta*TEMP + CPcor_SBE*PRES) / (1 + delta*TEMP_ADJUSTED + CPcor_new*PRES_ADJUSTED)',n,1); % First new entry
    %scientific_calib_coefficient.PSAL	= repmat('CPcor_new = xxx; CPcor_SBE = –9.57e–8; delta = 3.25e–6',n,1); % First new entry
    scientific_calib_coefficient.PSAL	= repmat(['CPcor_new = ',sprintf('%0.3g',CPcor_new),'; CPcor_SBE = ',sprintf('%0.3g',CPcor_SBE),'; delta = ',sprintf('%0.3g',d),'.'],n,1); % First new entry
    scientific_calib_comment.PSAL	= repmat('New conductivity computed by using a different CPcor value from that provided by Sea-Bird.',n,1); % First new entry
    scientific_calib.date.PSAL		= repmat(replace(datestr(now,30),'T',''),n,1); % First new entry

    % Create text-section for the report:
    fprintf(fid,'%s\n','\newpage');
    fprintf(fid,'%s\n','\subsection{Correcting deep ARVOR float pressure dependent conductivity bias}');
    fprintf(fid,'%s\n','For ARVOR-D floats, conductivity has to be corrected due to a too small conductivity cell compressibility term (\emph{CPcor}) in the manufacturer calibration of the SBE CTDs. Figure~\ref{prescondbias} shows the uncorrected and corrected salinities from the first cycle. ');
    %fprintf(fid,'%s\n','This re-computation salinity should be done as the first step in D-mode. However, to ensure the correct cycle is used for this calibration, it is done after the delayed-mode procedures for coordinates (Section~\ref{sec:DM-coordinates}). ');
    fprintf(fid,'%s\n','\begin{figure}[H]');
    fprintf(fid,'%s\n','\centering');
    fprintf(fid,'%s\n','\includegraphics[width=0.5\textwidth,natwidth=400,natheight=650]{\floatsource\WMOnum_prescondbias}');
    fprintf(fid,'%s%s%s%s%s\n','\caption{Float \WMOnum. Salinity profiles from the first cycle using the manufacturer calibration (i.e., raw data; $CPcor_{SBE} =$ ',num2str(CPcor_SBE),'~$dbar^{-1}$; blue) and the recommended compressibility term for sensor model ',ctdmodel,' ($CPcor_{new} =$ ',num2str(CPcor_new),'~$dbar^{-1}$; red).}');
    fprintf(fid,'%s\n','\label{prescondbias}');
    fprintf(fid,'%s\n','\end{figure}');
  else
    % Assuming that scientific_calib should be 'none' throughout:
    scientific_calib.equation.PSAL	= repmat('none',n,1); % First new entry 
    scientific_calib.coefficient.PSAL	= repmat('none',n,1); % First new entry
    scientific_calib.comment.PSAL	= repmat('none',n,1); % First new entry
    scientific_calib.date.PSAL  	= repmat(replace(datestr(now,30),'T',''),n,1); % First new entry
  end
  fclose(fid);
  
  % [] Check if this needs to be done earlier or later!

  
  
  
  
  %%%%%%% FINAL PHASE: Final plots, snippets, save for OWC, and alter R-files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %DENS = sw_dens(SAL,TEMP,PRES);	% Update density for the final plot:
  DENS = sw_pden(SAL,TEMP,PRES,0);	% Potential density instead

  % Finally plot the clean plot.
  plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'.eps']);
  
  % Create this new variable:
  PTMP = sw_ptmp(SAL,TEMP,PRES,0);

  %%%%%% Add flags to the Hov-Möller diagrams %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  figure(1001); 
  caxis(mima(PTMP));  % Set colorscale to the good data range
  find(TEMPqco=='4'); ho=line(X(ans),Y(ans)); set(ho,'marker','s','markersize',8,'color','w','linestyle','none');
  find(TEMPqcn==1);   hr=line(X(ans),Y(ans)); set(hr,'marker','.','markersize',8,'color','w','linestyle','none');
  find(TEMPqcn==4);   hn=line(X(ans),Y(ans)); set(hn,'marker','o','markersize',8,'color','w','linestyle','none');
  find(PRESqco=='4'); ho=line(X(ans),Y(ans)); set(ho,'marker','s','markersize',8,'color',[.5 .5 .5],'linestyle','none');
  find(PRESqcn==1);   hr=line(X(ans),Y(ans)); set(hr,'marker','.','markersize',8,'color',[.5 .5 .5],'linestyle','none');
  find(PRESqcn==4);   hn=line(X(ans),Y(ans)); set(hn,'marker','o','markersize',8,'color',[.5 .5 .5],'linestyle','none');
  print(gcf,'-dpng',[outfiles{I}(1:end-4),'_HMtemp.png']);
  figure(1002); 
  caxis(mima(SAL));
  find(SALqco=='4');  ho=line(X(ans),Y(ans)); set(ho,'marker','s','markersize',8,'color','w','linestyle','none');
  find(SALqcn==1);    hr=line(X(ans),Y(ans)); set(hr,'marker','.','markersize',8,'color','w','linestyle','none');
  find(SALqcn==4);    hn=line(X(ans),Y(ans)); set(hn,'marker','o','markersize',8,'color','w','linestyle','none');
  find(PRESqco=='4'); ho=line(X(ans),Y(ans)); set(ho,'marker','s','markersize',8,'color',[.5 .5 .5],'linestyle','none');
  find(PRESqcn==1);   hr=line(X(ans),Y(ans)); set(hr,'marker','.','markersize',8,'color',[.5 .5 .5],'linestyle','none');
  find(PRESqcn==4);   hn=line(X(ans),Y(ans)); set(hn,'marker','o','markersize',8,'color',[.5 .5 .5],'linestyle','none');
  print(gcf,'-dpng',[outfiles{I}(1:end-4),'_HMsal.png']);

  % %%%%%% Save for OWC and WRITE_D: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Add 360 to match the reference data:
  LONG<0; LONG(ans)=LONG(ans)+360; 
  %whos LAT LONG DATES PRES SAL TEMP PTMP PROFILE_NO CYCLE_NUMBER 
  save(outfiles{I}','LAT','LONG','DATES','PRES','SAL','TEMP','PTMP','PROFILE_NO','PROFILE_*_N','CYCLE_NUMBER','*qco','*qcn','Rfiles','scientific_calib');

  
  %%%%%%% The table of flags %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Go to that float's working dir and write snippets for qc-table:
  cd([my_working_dir,'DMQC',filesep,float_names{I}]);
  if ~isempty(comments)
    comments(end-1:end)='. '; comments=['New flags are due to: ',comments];  
  else
    comments='No new flags were found necessary.';
  end
  snippet(comments);
  fid=fopen('flagtable.tex','w');
  fprintf(fid,'%s\n','\begin{tabular}{|l|c|c|c|p{7cm}|}');
  fprintf(fid,'%s\n','\hline'); 
  fprintf(fid,'%s \\\\ \n',['Variable & RTQC flags (''4'') & reversed flags (''1'') & new flags (''4'') & Affected cycles (cycle numbers)']);
  fprintf(fid,'%s\n','\hline');
  fprintf(fid,'%s & %u & %u & %u & %s \\\\ \n','POS ',sum(POSqco=='4','all') ,sum(POSqcn==1,'all'), sum(POSqcn==4,'all'), snippet(CYCLE_NUMBER(find(any(POSqco=='4'|POSqcn==1|POSqcn==4)))));
  fprintf(fid,'%s & %u & %u & %u & %s \\\\ \n','JULD',sum(JULDqco=='4','all'),sum(JULDqcn==1,'all'),sum(JULDqcn==4,'all'),snippet(CYCLE_NUMBER(find(any(JULDqco=='4'|JULDqcn==1|JULDqcn==4)))));
  fprintf(fid,'%s & %u & %u & %u & %s \\\\ \n','PRES',sum(PRESqco=='4','all'),sum(PRESqcn==1,'all'),sum(PRESqcn==4,'all'),snippet(CYCLE_NUMBER(find(any(PRESqco=='4'|PRESqcn==1|PRESqcn==4)))));
  fprintf(fid,'%s & %u & %u & %u & %s \\\\ \n','PSAL',sum(SALqco=='4','all') ,sum(SALqcn==1,'all') ,sum(SALqcn==4,'all'), snippet(CYCLE_NUMBER(find(any(SALqco=='4'|SALqcn==1|SALqcn==4)))));
  fprintf(fid,'%s & %u & %u & %u & %s \\\\ \n','TEMP',sum(TEMPqco=='4','all'),sum(TEMPqcn==1,'all'),sum(TEMPqcn==4,'all'),snippet(CYCLE_NUMBER(find(any(TEMPqco=='4'|TEMPqcn==1|TEMPqcn==4)))));
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

