% The subdirectory /data is organised as follows:
%
% /data/float_source/
% contains .mat files with the original data from the floats.
% /data/float_mapped/
% contains .mat files with the objective estimates at the float profile locations and observed θ levels.
% /data/float_calib/
% contains .mat files with the calibration constants and error estimates.
% /data/float_plots/
% contains diagnostic plots from the calibration.
% /data/constants/
% contains coastdat.mat, wmo_boxes.mat, and TypicalProfileAroundSAF.mat.
% /data/climatology/historical_ctd, /historical_bot, /argo_profiles
% are where you put your reference data. Reference data are saved as .mat files in 10×10 degree WMO boxes. Please refer to README_prepare_refdbase_ow.pdf for their data format.
 
% ----- Download float data manually: ------------------
%
% cd ~/Downloads/ARGO/6902550/profiles
% lftp ftp.ifremer.fr/ifremer/argo/dac/coriolis/6902550/profiles
									% Only D and R files???
%
% ------------------------------------------------------

% Paths and filenames:
clear all
float_names={'6902550','5904988'};			% Necessary for ow_calibration.m. Do not change!
float_names={'6902550'};
float_names={'5904988'};
%float_names={'3902103'};
% Root folders:
down_dir=strcat('~/Downloads/ARGO/',float_names,filesep)
rootdirin=strcat(down_dir,'Rfiles/')
rootdirout=strcat(down_dir,'Dfiles/')
% Matlabow folders:
source_dir='~/Arkiv/data/matlabow/float_source/';
float_dirs=repmat({'project_NorARGO/'},1,length(float_names)); % Necessary for ow_calibration.m. Do not change!
% Files:
infiles=strcat(down_dir,filesep,float_names,{'_prof.nc'}) 	% Using _prof ncfiles.
outfiles=strcat(source_dir,float_dirs,float_names,'.mat')

% infiles=strcat(rootdirin,filesep,float_names,{'.nc'});	% For using profile-files as infiles for matlabow (maybe later).
% for I=1:length(float_names) 
  % %infiles=edir(rootdirin{I},'nc',0,0);  
  % files=dir(rootdirin{I});
  % char(files.name);find(ans(:,1)=='D' | ans(:,1)=='R'); char(files(ans).name);
  % %ans(:,1)='D';
  % filenames=sortrows(ans)
  % mkdir(rootdirout{I});
  % system(['cp ',rootdirin{I},filesep,'D',float_names{I},'_*.nc ',rootdirout{I},filesep])
  % system(['cp ',rootdirin{I},filesep,'R',float_names{I},'_*.nc ',rootdirout{I},filesep])
  % return

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
  j=find(POSqc==4 | JULDqc==4 & -90>LAT & LAT>90);  
  LONG(j)=NaN; LAT(j)=NaN; DATES(j)=NaN; PRES(:,j)=NaN; SAL(:,j)=NaN; TEMP(:,j)=NaN;

  % QC flag remove measurements:
  PRESqc=ncread(infiles{I},'PRES_QC'); PRES(PRESqc==4)=NaN;
  SALqc=ncread(infiles{I},'PSAL_QC'); SAL(SALqc==4)=NaN;
  TEMPqc=ncread(infiles{I},'TEMP_QC'); TEMP(TEMPqc==4)=NaN;
  
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
  
  % Monotonically increasing pressure only:
  [logical(zeros(1,n));diff(PRES,1,1)<=0]; 
  if any(ans,'all')
    %PRES(ans)=NaN; SAL(ans)=NaN; TEMP(ans)=NaN;
    warning(['Non-monotonic pressure in profile(s) ',int2str(find(any(ans))),' !']);	% Only warning
  end
  
  % Check for density inversions:
  DENS = sw_dens(SAL,TEMP,PRES);
  [logical(zeros(1,n));diff(DENS,1,1)<=0];
  if any(ans,'all')
    %PRES(ans)=NaN; SAL(ans)=NaN; TEMP(ans)=NaN; 
    warning(['Density inversions in profile(s) ',int2str(find(any(ans))),' !']);	% Only warning
  end
  % Upwards and downwards?
  return

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

