% This script loads float nc-files from your download folder, checks
% and prepares float data for ow_calibration. Downloading from server
% can also be done here.

clear all
init_dmqc; % Paths and filenames, and ftp download.

if logical(0) % ----- SWITCH ON IN ORDER TO DOWNLOAD NEW FLOAT FILES AND R-FILES ------
  for I=1:length(download_dir)
    mkdir(download_dir{I}); 
    system(['lftp -e ''lcd ',download_dir{I},' ; get ',float_names{I},'_prof.nc ; bye'' ftp.ifremer.fr/ifremer/argo/dac/coriolis/',float_names{I},'/']);
    mkdir(rootdirin{I}); 
    system(['lftp -e ''lcd ',rootdirin{I},' ; mget R',float_names{I},'_*.nc ; bye'' ftp.ifremer.fr/ifremer/argo/dac/coriolis/',float_names{I},'/profiles/']);
  end
end % --------------------------------------------------------------------------------
% R-files must be downloaded at the same time, or else there might be
% mismatch between number of cal parameters and R-files.

% INGEST FLOAT FILES INTO THE OWC SYSTEM:

for I=1:length(float_names)
  
  % Reading and some massage:
  LONG=ncread(infiles{I},'LONGITUDE')'; 
  LONG<0; LONG(ans)=LONG(ans)+360; % Add 360 to match the reference data.
  LAT=ncread(infiles{I},'LATITUDE')';
  JULD=ncread(infiles{I},'JULD'); REFERENCE_DATE_TIME=ncread(infiles{I},'REFERENCE_DATE_TIME');
  REFERENCE_DATE_TIME'; DATES=dyear(datenum([ans(1:8),'T',ans(9:14)],'yyyymmddTHHMMSS')+JULD');
  PRES=ncread(infiles{I},'PRES');
  SAL=ncread(infiles{I},'PSAL');
  TEMP=ncread(infiles{I},'TEMP');
  CYCLE_NUMBER=ncread(infiles{I},'CYCLE_NUMBER')';
  DIRECTION=ncread(infiles{I},'DIRECTION')';

  % QC flag exclude rows and columns:
  POSqc=ncread(infiles{I},'POSITION_QC')'; JULDqc=ncread(infiles{I},'JULD_QC')'; 
  %j=find(POSqc~=4 & JULDqc~=4 & -90<=LAT & LAT<=90);  
  %LONG=LONG(j); LAT=LAT(j); DATES=DATES(j); CYCLE_NUMBER=CYCLE_NUMBER(j); PRES=PRES(:,j); SAL=SAL(:,j); TEMP=TEMP(:,j);
  % NaN out instead:
  j=find(POSqc=='4' | JULDqc==4 & -90>LAT & LAT>90);  
  LONG(j)=NaN; LAT(j)=NaN; DATES(j)=NaN; PRES(:,j)=NaN; SAL(:,j)=NaN; TEMP(:,j)=NaN;

  % QC flag remove measurements:
  PRESqc=ncread(infiles{I},'PRES_QC'); PRES(PRESqc=='4')=NaN;
  SALqc=ncread(infiles{I},'PSAL_QC');   SAL( SALqc=='4')=NaN;
  TEMPqc=ncread(infiles{I},'TEMP_QC'); TEMP(TEMPqc=='4')=NaN;

  % Sort cycle number:
  [ans,IA]=sort(CYCLE_NUMBER); % C = A(IA);
  if ~all(ans==CYCLE_NUMBER)
    LONG=LONG(IA); LAT=LAT(IA); DATES=DATES(IA); CYCLE_NUMBER=CYCLE_NUMBER(IA); PRES=PRES(:,IA); SAL=SAL(:,IA); TEMP=TEMP(:,IA);
    warning('Cycle numbers not in succession! (sorted)');
    %CYCLE_NUMBER=ans;
  end
  [m,n]=size(PRES);
  PROFILE_NO=1:n;
  % Maybe consider using DIRECTION in the sorting too, so it is
  % always, e.g., DADAAA...
  
  DENS = sw_dens(SAL,TEMP,PRES); % For the plotting of warning plots.
  zerow=zeros(1,n);			% Row of zeros

  % Pressure increasing test / monotonically increasing pressure test:
  logical([zerow;diff(PRES,1,1)<=0]); 
  ans(PRES<=MAP_P_EXCLUDE)=logical(0); % Do not be concerned with monotonicity above.
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
    clear pnb jnb nbt
  else
    system(['rm -f ',outfiles{I}(1:end-4),'_pressure_warning.eps']);
  end

  
  % Spike tests (RTQC double check and stricter DMQC check for some):
  % Test value = | V2 – (V3 + V1)/2 | – | (V3 – V1) / 2 |
  % according to EuroGOOS,  where V2 is the measurement being tested
  % as a spike, and V1 and V3 are the values above and below. 
  testvalue = [zerow; abs(SAL(2:end-1,:)-(SAL(3:end,:)+SAL(1:end-2,:))/2) - abs((SAL(3:end,:)-SAL(1:end-2,:))/2) ; zerow];
  %[PRES<500 & testvalue>0.9 |  PRES>=500 & testvalue>0.3]; % The RTQC9 limits 
  [PRES<500 & testvalue>0.9 |  PRES>=500 & testvalue>0.02]; % By experience clear spikes in the Nordic Seas
  if any(ans,'all')
    jnb=find(any(ans)); snb=ans; nbt='Salinity spike(s)';
    warning([nbt,' in float ',float_names{I},' profile(s): ',int2str(PROFILE_NO(jnb)),' ! (removed)']);
    plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_S-spike_warning.eps']);
    SAL(snb)=NaN; % REMEMBER TO ALSO FLAG IT
    clear snb jnb nbt
  else
    system(['rm -f ',outfiles{I}(1:end-4),'_S-spike_warning.eps']);
  end
  testvalue = [zerow; abs(TEMP(2:end-1,:)-(TEMP(3:end,:)+TEMP(1:end-2,:))/2) - abs((TEMP(3:end,:)-TEMP(1:end-2,:))/2) ; zerow];
  [PRES<500 & testvalue>6 |  PRES>=500 & testvalue>2];
  if any(ans,'all')
    jnb=find(any(ans)); tnb=ans; nbt='Temperature spike(s)';
    warning([nbt,' in float ',float_names{I},' profile(s): ',int2str(PROFILE_NO(jnb)),' ! (removed)']);
    plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_T-spike_warning.eps']);
    TEMP(tnb)=NaN; % REMEMBER TO ALSO FLAG IT
    clear tnb jnb nbt
  else
    system(['rm -f ',outfiles{I}(1:end-4),'_T-spike_warning.eps']);
  end


  % Gradient test
  % Test value = | V2 − (V3 + V1)/2 |
  % where V2 is the measurement being tested, and V1 and V3 are the values above and below.
  testvalue = [zerow ; abs(SAL(2:end-1,:)-(SAL(3:end,:)+SAL(1:end-2,:))/2) ; zerow];
  [PRES<500 & testvalue>1.5 |  PRES>=500 & testvalue>0.5]; % The RTQC9 limits 
  %[PRES<500 & testvalue>1.5 |  PRES>=500 & testvalue>0.5]; % By experience in the Nordic Seas  
    if any(ans,'all')
    jnb=find(any(ans));  snb=ans; nbt='Salinity gradient(s)';
    warning([nbt,' in float ',float_names{I},' profile(s): ',int2str(PROFILE_NO(jnb)),' ! (removed)']);
    plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_S-gradient_warning.eps']);
    SAL(snb)=NaN; % REMEMBER TO ALSO FLAG IT
    clear snb jnb nbt
  else
    system(['rm -f ',outfiles{I}(1:end-4),'_S-gradient_warning.eps']);
  end
  testvalue = [zerow ; abs(TEMP(2:end-1,:)-(TEMP(3:end,:)+TEMP(1:end-2,:))/2) ; zerow];
  [PRES<500 & testvalue>9 |  PRES>=500 & testvalue>3]; % The RTQC9 limits 
  if any(ans,'all')
    jnb=find(any(ans)); tnb=ans; nbt='Temperature gradient(s)';
    warning([nbt,' in float ',float_names{I},' profile(s): ',int2str(PROFILE_NO(jnb)),' ! (removed)']);
    plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_T-gradient_warning.eps']);
    TEMP(tnb)=NaN; % REMEMBER TO ALSO FLAG IT
    clear tnb jnb nbt 
  else
    system(['rm -f ',outfiles{I}(1:end-4),'_T-gradient_warning.eps']);
  end

  % Density inversion test:
  pres=PRES; pres(isnan(SAL)|isnan(TEMP))=NaN;
  PR = pres - [zerow; nandiff(pres)/2];	% Reference pressure
  DENS = sw_pden(SAL,TEMP,pres,PR);	% Potential density
  downw=nandiff(DENS)<-0.03;		% Downward test
  pres([downw;logical(zerow)])=NaN;	% Remove found
  PR = pres - [zerow; nandiff(pres)/2];	% Reference pressure again
  dens = sw_pden(SAL,TEMP,pres,PR);	% Potential density again
  upw=nandiff(dens)<-0.03;			% Upward test
  logical([downw;zerow]+[zerow;upw]);	% Logical for bad data
  if any(ans,'all')
    jnb=find(any(ans)); inb=ans; nbt='Density inversions';
    warning([nbt,' in float ',float_names{I},' profile(s): ',int2str(PROFILE_NO(jnb)),' !']);	% Only warning
    plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_inversion_warning.eps'])
    SAL(inb)=NaN; TEMP(inb)=NaN; DENS(inb)=NaN;  % REMEMBER TO ALSO FLAG S&T
    clear inb jnb nbt
  else
    system(['rm -f ',outfiles{I}(1:end-4),'_inversion_warning.eps']);
  end

  % Finally plot the clean plot.
  plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'.eps']);

  
  % Create this new variable:
  PTMP = sw_ptmp(SAL,TEMP,PRES,0);
  
  %whos LAT LONG DATES PRES SAL TEMP PTMP PROFILE_NO CYCLE_NUMBER 
  save(outfiles{I}','LAT','LONG','DATES','PRES','SAL','TEMP','PTMP','PROFILE_NO','CYCLE_NUMBER');
end


% Checklist:
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

