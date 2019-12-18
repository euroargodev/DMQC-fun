% This script loads float nc-files from your download folder, checks
% and prepares float data for ow_calibration.

clear all
init_dmqc; % Paths and filenames, and ftp download.

% ----- Download float data manually: ------------------
% OBSOLETE. This is now done in init_dmqc.
% cd ~/Downloads/ARGO/
% lftp ftp.ifremer.fr/ifremer/argo/dac/coriolis/
% ------------------------------------------------------

for I=1:length(float_names)
  
  % Reading and some massage:
  LONG=ncread(infiles{I},'LONGITUDE')'; LONG<0; LONG(ans)=LONG(ans)+360;
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
    CYCLE_NUMBER=ans;
  end
  [m,n]=size(PRES);
  PROFILE_NO=1:n;
  % Maybe consider using DIRECTION in the sorting too, so it is
  % always, e.g., DADAAA...
  
  DENS = sw_dens(SAL,TEMP,PRES);

  % Spike tests (RTQC double check and stricter DMQC check for some):
  % Test value = | V2 – (V3 + V1)/2 | – | (V3 – V1) / 2 |
  % according to EuroGOOS,  where V2 is the measurement being tested
  % as a spike, and V1 and V3 are the values above and below. 
  testvalue = [zeros(1,n); abs(SAL(2:end-1,:)-(SAL(3:end,:)+SAL(1:end-2,:))/2) - abs((SAL(3:end,:)-SAL(1:end-2,:))/2) ; zeros(1,n)];
  %[PRES<500 & testvalue>0.9 |  PRES>=500 & testvalue>0.3]; % The RTQC9 limits 
  [PRES<500 & testvalue>0.9 |  PRES>=500 & testvalue>0.02]; % By experience clear spikes in the Nordic Seas
  if any(ans,'all')
    jnb=find(any(ans)); nb=ans; nbt='Salinity spike(s)';
   find(nb)
    warning([nbt,' in float ',float_names{I},' profile(s): ',int2str(jnb),' ! (removed)']);
    plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_S-spike_warning.eps'])
    SAL(nb)=NaN; % REMEMBER TO ALSO FLAG IT
  else
    clear nb
    system(['rm -f ',outfiles{I}(1:end-4),'_S-spike_warning.eps']);
  end
  testvalue = [zeros(1,n); abs(TEMP(2:end-1,:)-(TEMP(3:end,:)+TEMP(1:end-2,:))/2) - abs((TEMP(3:end,:)-TEMP(1:end-2,:))/2) ; zeros(1,n)];
  [PRES<500 & testvalue>6 |  PRES>=500 & testvalue>2];
  if any(ans,'all')
    jnb=find(any(ans)); nb=ans; nbt='Temperature spike(s)';
    warning([nbt,' in float ',float_names{I},' profile(s): ',int2str(jnb),' ! (removed)']);
    plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_T-spike_warning.eps']);
    TEMP(nb)=NaN; % REMEMBER TO ALSO FLAG IT
  else
    clear nb
    system(['rm -f ',outfiles{I}(1:end-4),'_T-spike_warning.eps']);
  end

  DENS = sw_dens(SAL,TEMP,PRES);

  % Monotonically increasing pressure only:
  [logical(zeros(1,n));diff(PRES,1,1)<=0]; 
  ans(PRES<=MAP_P_EXCLUDE)=logical(0); % Do not be concerned with monotonicity above.
  if any(ans,'all')
    %PRES(ans)=NaN; SAL(ans)=NaN; TEMP(ans)=NaN;
    jnb=find(any(ans)); nb=ans; nbt='Non monotonic pressure';
    warning(['Non-monotonic pressure in float ',float_names{I},' profile(s): ',int2str(jnb),' !']);	% Only warning
    plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_pressure_warning.eps'])
  else
    system(['rm -f ',outfiles{I}(1:end-4),'_pressure_warning.eps']);
  end
  
  % Check for density inversions:
  [logical(zeros(1,n));diff(DENS,1,1)<=0]; 
  ans(PRES<=MAP_P_EXCLUDE)=logical(0); % Do not be concerned with inversions above.
  if any(ans,'all')
    %PRES(ans)=NaN; SAL(ans)=NaN; TEMP(ans)=NaN; 
    jnb=find(any(ans)); nb=ans; nbt='Inversions';
    warning(['Density inversions in float ',float_names{I},' profile(s): ',int2str(jnb),' !']);	% Only warning
    plot_profiles; print(gcf,'-depsc',[outfiles{I}(1:end-4),'_inversion_warning.eps'])
  else
    system(['rm -f ',outfiles{I}(1:end-4),'_inversion_warning.eps']);
  end
  % Upwards and downwards?

  % Create this new variable:
  PTMP = sw_ptmp(SAL,TEMP,PRES,0);
  
  whos LAT LONG DATES PRES SAL TEMP PTMP PROFILE_NO CYCLE_NUMBER
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

