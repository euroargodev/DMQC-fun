% WRITE_D script writes the results of DMQC to D-files. 
%
% DMQC-fun v0.9.
% J. Even Ø. Nilsen, Ingrid M. Angel-Benavides, Birgit Klein, Malgorzata Merchel, and Kjell Arne Mork.
% Last updated: Wed May 24 13:40:44 2023 by jan.even.oeie.nilsen@hi.no

% This is the last thing you do in DMQC.
%
% You set which floats to operate on etc. in INIT_DMQC. 
%
% No editing in this file!
%
% ------------------------------------------------------------------------------
%
% Some notes: The syncronisity between columns and files is secured by
% PREPARE_FLOATS and WRITE_D both being based on the accessible files,
% hence the columns of the matrices saved by PREPARE_FLOATS corresponds
% to the files. Furthermore, the fil list is the same, as it is the list
% made by PREPARE_FLOATS that is used by WRITE_D. Data for ascending and
% descending profiles are saved in separate directories, and run
% separately. RUN_OW_CALIBRATION only operates on ascending profiles,
% but the results, ... , are stored in ...
%
% Herein some passages from the QC manual is quoted. These are
% commented with %W
%
% ------------------------------------------------------------------------------

clear all; close all
init_dmqc; % Paths and filenames.

% For write_D, we always run both directions at the same time:
direction={'A','D'};

for II=1:length(direction)	% Loop direction

  for I=1:length(download_dir)	% Loop floats

    load(outfiles{I}); % data from PREPARE_FLOATS
	% Listed only what is used herein     
	% 1: 'PRES','TEMP','PTMP',...
	% 2: 'JULD','CYCLE_NUMBER',...
	% 3: 'PROFILE_*_N','Rfiles','scientific_calib','PSAL_ADJUSTED_Cnew',...
	% 4: 'PSAL'); 
	% First row is for OWC, as well as OPERATOR_CPCOR_NEW and WRITE_D;
	% Second row only CYCLE_NUMBER is used in write_D;
	% Third row and PSAL is for WRITE_D;
	% Fourth row is for OPERATOR_CPCOR_NEW.
	% Needed for both ascending and descending case.
	
    load([calib_dir,'calseries_',float_names{I}],'calseries'); % settings for OWC
	% breaks                0x0                
	% calib_profile_no      1x113              
	% calseries*            1x113              
	% max_breaks            1x1                
	% use_percent_gt        1x1                
	% use_pres_gt           0x0                
	% use_pres_lt           0x0                
	% use_theta_gt          0x0                
	% use_theta_lt          0x0                
	% *: the only objects used in this script.
	% Needed for both ascendiing and descending case.
	
    load([calib_dir,'cal_',float_names{I}]);	% results of OWC
	% PROFILE_NO               1x113              
	% cal_COND*             2500x113              
	% cal_COND_err*         2500x113              
	% cal_SAL*              2500x113  [√]            
	% cal_SAL_err*          2500x113  [√]            
	% fbreaks                  0x0                
	% fcoef                    1x1                
	% pcond_factor*            1x113              
	% pcond_factor_err*        1x113              
	% sta_SAL*              2500x113              
	% sta_SAL_err*          2500x113              
	% sta_mean                 1x113              
	% sta_rms                  1x113              
	% time_deriv               1x113              
	% time_deriv_err           1x113            
	% *: the objects used in this script.
	% Needed for both ascendiing and descending case.
	

    % ---- Assign float_working_dir (and calibration of descending profiels)  ----
    switch direction{II}
     case 'A',
      disp(['Writing ascending profiles from float ',float_names{I},'.']);
      float_working_dir=[my_working_dir,'DMQC',filesep,float_names{I}];
     case 'D'
      disp(['Writing DESCENDING profiles from float ',float_names{I},'.']);
      float_working_dir=[my_working_dir,'DMQC',filesep,float_names{I},filesep,'descending'];
      if ~exist(float_working_dir,'dir')
	error('Descending profiles not processed yet!')
      end
      % Keep the cycle numbers for ascending profiles:    %load(outfiles{I},'CYCLE_NUMBER'); 
      CYCLE_NUMBER_a=CYCLE_NUMBER;
      if any(diff(CYCLE_NUMBER_a)<1), error('CYCLE_NUMBER_a from ',outfiles{I},' not increasing!'); end
      % Redirect to load data from descending outfile:
      outfiles{I}=[outfiles{I}(1:end-4),'D',outfiles{I}(end-3:end)];
      if exist(outfiles{I},'file')
	load(outfiles{I}); 
      else
	['There are no direction ''',direction{II},''' files for float ',float_names{I},'. '];
	disp(ans);
	continue 
      end
      if any(diff(CYCLE_NUMBER)<1), error('CYCLE_NUMBER from ',outfiles{I},' not increasing!'); end
      % CALIBRATION of descending profiles is done by applying the
      % corrections for the ascending profile in the same cycle. Hence, we
      % match the cycle numbers from ascending and descending, and map the
      % descending ones to the ascending:
      calseries_a=calseries;				% The loaded calseries is for ascending profiles
      calseries=zeros(size(CYCLE_NUMBER));		% New calseries for the (current) descending profiles
      pcond_factor_a=pcond_factor;			% The loaded pcond_factor is for ascending profiles
      pcond_factor=ones(size(CYCLE_NUMBER));		% New pcond_factor for the (current) descending profiles
      pcond_factor_err_a=pcond_factor_err;		% The loaded pcond_factor_err is for ascending profiles
      pcond_factor_err=zeros(size(CYCLE_NUMBER));	% New pcond_factor_err for the (current) descending profiles
      [~,ia,ib] = intersect(CYCLE_NUMBER,CYCLE_NUMBER_a); % Indices for matching cycle numbers   
      calseries(ia)=calseries_a(ib);			% Assign the corresponding calseries numbers 
      cal_action{I}=cal_action{I}(unique(calseries(calseries>0))); % Pick only the relevant cal_actions
      pcond_factor(ia)=pcond_factor_a(ib);		% Assign the corresponding pcond_factor numbers 
      pcond_factor_err(ia)=pcond_factor_err_a(ib);	% Assign the corresponding pcond_factor_err numbers 
      clear *_a ia ib
      % Calculate descending cal_SAL and cal_SAL_err to be used below, the
      % same way as OWC's CALCULATE_PIECEWISEFIT:
      COND = sw_c3515*sw_cndr(SAL,PTMP,0);
      cal_COND = pcond_factor.*COND;
      cal_SAL = sw_salt(cal_COND/sw_c3515,PTMP,0);
      cal_COND_err = pcond_factor_err.*COND;
      cal_SAL_err = abs(cal_SAL-sw_salt((cal_COND+cal_COND_err)/sw_c3515,PTMP,0));
    end
    cd(float_working_dir);

  
    % ---- Create D-files ----
    % Rfiles object have been loaded from PREPARE_FLOATS output above.
    % We write both ascending and descening profiles in this script.
    % Check for unknown file types:
    if ~all(contains(Rfiles,{'/R','/D'}))		
      error('Other files than R and D files exists! Maybe A files? Must recode!');
    end
    % Clear out existing directory of D-files:
    if  ~exist(rootdirout{I},'dir'), mkdir(rootdirout{I});  
    elseif II==1, rmdir(rootdirout{I},'s'); mkdir(rootdirout{I});  
    end
    % Rename Rfiles to Dfiles, change path, and copy from one to other:
    Dfiles=replace(replace(Rfiles,'//R','//D'),rootdirin{I},rootdirout{I});
    command=strcat({'cp '},Rfiles,{' '},Dfiles);
    for i=1:length(command), system(command{i}); end
    % Count files and test vs. matrix data:
    FN=length(Dfiles);
    if FN~=size(PRES,2), error('Number of files does not match number of columns in data!'); end
    
    % ---- Load the flags from DMQC ----
    load(['savedflags_',float_names{I}],'*qcn','*qco');
	% JULDqcn        1x17               
	% JULDqco        1x17               
	% POSqcn         1x17               
	% POSqco         1x17               
	% PRESqcn      700x17               
	% PRESqco      700x17               
	% PSALqcn      700x17               
	% PSALqco      700x17               
	% TEMPqcn      700x17               
	% TEMPqco      700x17               
  
    % ---- Prepare the data and flag matrices ----
    m=size(PRES,1)-1; % Remove the row added to make sure every profile ends with a NaN.
    PRES                = PRES(1:m,:);                
    PRESqcn             = PRESqcn(1:m,:);             
    PRESqco             = PRESqco(1:m,:);             
    PSAL                = PSAL(1:m,:);                
    PSALqcn             = PSALqcn(1:m,:);             
    PSALqco             = PSALqco(1:m,:);             
    PTMP                = PTMP(1:m,:);                
    SAL                 = SAL(1:m,:);                 
    TEMP                = TEMP(1:m,:);                
    TEMPqcn             = TEMPqcn(1:m,:);             
    TEMPqco             = TEMPqco(1:m,:);             
    cal_COND            = cal_COND(1:m,:);            
    cal_COND_err        = cal_COND_err(1:m,:);        
    cal_SAL             = cal_SAL(1:m,:);             
    cal_SAL_err         = cal_SAL_err(1:m,:);         
    sta_SAL             = sta_SAL(1:m,:);             
    sta_SAL_err         = sta_SAL_err(1:m,:);
    if ~isempty(PSAL_ADJUSTED_Cnew), PSAL_ADJUSTED_Cnew  = PSAL_ADJUSTED_Cnew(1:m,:); end
    % whos PRES* TEMP* PSAL* SAL PTMP sta_SAL* cal_SAL* cal_COND*
    %
    % Initiate variables for NetCDF, back from the OWC renaming in PREPARE_FLOATS (PSAL is renamed further below):
    PRES_ADJUSTED=PRES;		TEMP_ADJUSTED=TEMP;				% Carry over the old data
    PRES_ADJUSTED_QC=PRESqco;	TEMP_ADJUSTED_QC=TEMPqco;			% Carry over the old flags
    PRES_QC=PRESqco;              TEMP_QC=TEMPqco;	PSAL_QC=PSALqco;	% Rename all old flags to fit netcdf names
    % Delete obsolete variables:
    clear *qco		% Obsolete from here on
    clear TEMP		% Obsolete from here on
    % PRES		% The original PRES is used for PRES_ADJUSTED_ERROR below
    clear PTMP SAL	% Obsolete after calculating calibration for descending profiles (above) 
    % *qcn		% DMQC flags are assigned to <PARAMETER>_ADJUSTED_QC below
  
  
    % ---- Add the results of OWC to objects ----
    [PSAL_ADJUSTED,PSAL_ADJUSTED_ERROR]=deal(nan(size(PSAL)));
    % For two cases below old flags are just carried over, while for the
    % cases where data are bad and unadjustable, the profiles will be
    % flagged '4'. So to carry over flags for profiles OWC ignores, we
    % duplicate old flags initially:
    PSAL_ADJUSTED_QC=PSAL_QC;	% RTQC flags added here (DMQC flags
    % are added further below). The calseries contains sequences of
    % numbers (e.g., 1 1 1 2 2 3 3 3) indicating groups of profiles the
    % calibration actions are to be applied to. Hence, we need to check
    % consistency with your cal_action vector (0 means no action; 1 means
    % adjust according to calibration; 4 means data are bad and
    % unadjustable):
    Na=length(cal_action{I});
    if length(unique(calseries(calseries>0)))~=Na, error('Number of specified actions differ from parts of calseries!'); end
    %
    for i=1:Na			% Loop cal_actions (i.e., my decisions)
      jj=find(calseries==i);							% Indices to profiles for i-th action 
      dPJ=contains(string(scientific_calib.equation.PRES(jj,:)),'PRES - dP');	% True at profiles of these with SP correction
      switch cal_action{I}(i)							% What is i-th action?
       case 0			% Data are good; no PSAL adjustment has been applied by OWC:
	PSAL_ADJUSTED_ERROR(:,jj) = max(cal_SAL_err(:,jj),0.01);	% Yes, also when no adjustment is done.[]
	if isempty(PSAL_ADJUSTED_Cnew)
	  PSAL_ADJUSTED(:,jj) = PSAL(:,jj); % PSAL maybe SP corrected, or original PSAL, with NaNs assigned by PREPARE_FLOATS
	  'PSAL_ADJUSTED = PSAL. ';					% No CPcorr, no SP, no drift
		scientific_calib.equation.PSAL(jj(~dPJ),end+[1:size(ans,2)]) = repmat(ans,length(jj(~dPJ)),1);
	  'PSAL_ADJUSTED = PSAL(COND, TEMP_ADJUSTED, PRES_ADJUSTED). ';	% No CPcorr, has SP, no drift'
		scientific_calib.equation.PSAL(jj(dPJ),end+[1:size(ans,2)]) = repmat(ans,length(jj(dPJ)),1);
	else
	  % PSAL_ADJUSTED_Cnew formula involving PRES_ADJUSTED etc. is found in manual.
	  PSAL_ADJUSTED(:,jj) = PSAL_ADJUSTED_Cnew(:,jj); 
	  'PSAL_ADJUSTED = PSAL_ADJUSTED_Cnew. ';				% Has CPcorr, regardless of SP, no drift
		scientific_calib.equation.PSAL(jj,end+[1:size(ans,2)]) = repmat(ans,length(jj),1);
	end
	% Coefficients and comments regarding OWC:
	'none';
		scientific_calib.coefficient.PSAL(jj,end+[1:size(ans,2)]) = repmat(ans,length(jj),1);
        'No significant salinity offset or drift detected. The quoted error is max[0.01, statistical uncertainty] in PSS-78. ';
		scientific_calib.comment.PSAL(jj,end+[1:size(ans,2)]) = repmat(ans,length(jj),1);
        %W PSAL_ADJUSTED = PSAL (original value);
	%W PSAL_ADJUSTED_QC = ‘1’ or ‘2’;
	%W PSAL_ADJUSTED_ERROR = maximum [statistical uncertainty, 0.01].  
       case 1			% Data show sensor drift or offset; adjustment has been applied:
	PSAL_ADJUSTED(:,jj) = cal_SAL(:,jj);	% OWC carries over the CPnew corrections together with its adjustment, in cal_SAL
						% But the scientifc cal equation becomes different.
        PSAL_ADJUSTED_ERROR(:,jj) = max(cal_SAL_err(:,jj),0.01);
	if isempty(PSAL_ADJUSTED_Cnew)
	  'PSAL_ADJUSTED = PSAL + dS. ';						% No CPcorr, no SP, has drift
      		scientific_calib.equation.PSAL(jj(~dPJ),end+[1:size(ans,2)]) = repmat(ans,length(jj(~dPJ)),1);
	  'PSAL_ADJUSTED = PSAL(COND, TEMP_ADJUSTED, PRES_ADJUSTED) + dS. ';	% No CPcorr, has SP, has drift
		scientific_calib.equation.PSAL(jj(dPJ),end+[1:size(ans,2)]) = repmat(ans,length(jj(dPJ)),1);
	else
	  % PSAL_ADJUSTED_Cnew formula involving PRES_ADJUSTED etc. is found in manual.
	  'PSAL_ADJUSTED = PSAL_ADJUSTED_Cnew + dS. ';				% Has CPcorr, regardless of SP, has drift
		scientific_calib.equation.PSAL(jj,end+[1:size(ans,2)]) = repmat(ans,length(jj),1);
	end
	% Coefficients and comments regarding OWC:
	% 'Vertically averaged dS = 0.012 +/- 0.003';
	strcat('Vertically averaged dS = ',num2str(mean([cal_SAL(:,jj)-PSAL(:,jj)],1,'omitnan')','%5.3f'),' +/- ',num2str(mean(cal_SAL_err(:,jj),1,'omitnan')','%5.3f'),'. ');
		scientific_calib.coefficient.PSAL(jj,end+[1:size(ans,2)]) = ans;
        'Significant salinity sensor offset or drift detected. OW least squares fit adopted. The quoted error is max[0.01, statistical uncertainty] in PSS-78. ';
		scientific_calib.comment.PSAL(jj,end+[1:size(ans,2)]) = repmat(ans,length(jj),1);
		scientific_calib.date.PSAL(jj,:) = repmat(replace(datestr(now,30),'T',''),length(jj),1); 
        %W PSAL_ADJUSTED = original value + adjustment recommended by statistical analysis, or adjustment provided by PI;
	%W PSAL_ADJUSTED_QC = ‘1’ or ‘2’;
	%W PSAL_ADJUSTED_ERROR = maximum [ (Σadjustment_error2)1/2, 0.01]
	%W where “adjustment_error” is the uncertainty from each type of adjustment applied to PSAL. These
	%W can be statistical uncertainty from sensor drift adjustment, uncertainty from conductivity cell
	%W thermal mass adjustment, etc.
       case 4			% Data are bad and unadjustable:
	PSAL_ADJUSTED(:,jj) = NaN; 
	PSAL_ADJUSTED_QC(:,jj) = '4';
	PSAL_ADJUSTED_ERROR(:,jj) = NaN;
	'none';	scientific_calib.equation.PSAL(jj,end+[1:size(ans,2)]) = repmat(ans,length(jj),1);
	'none';	scientific_calib.coefficient.PSAL(jj,end+[1:size(ans,2)]) = repmat(ans,length(jj),1);
	'Salinity data are bad and unadjustable. ';
		scientific_calib.comment.PSAL(jj,end+[1:size(ans,2)]) = repmat(ans,length(jj),1);
		scientific_calib.date.PSAL(jj,:) = repmat(replace(datestr(now,30),'T',''),length(jj),1); 
        %W "PSAL_ADJUSTED = FillValue;"			
	%W "PSAL_ADJUSTED_QC = ‘4’; PSAL_QC = '4';"	
	%W "PSAL_ADJUSTED_ERROR = FillValue."
	%W   scientific_calib.equation.PSAL = 'none';
	%W   scientific_calib.coefficient.PSAL = 'none';
	%W   scientific_calib.comment.PSAL = 'Salinity data are bad and unadjustable because of the Druck oil microleak problem.';
	% NaNs are replaced by FillValues below.
      end % What i-th action switch
    end % Loop cal_actions
    % This naive adding to the scientific calib info, often results in
    % superfluous 'none' in the texts. This will be cleaned up below (3.6.2).


    % Some times OWC removes data above a set of NaN values; flag these
    % as bad (qc='4'):
    PSALqcn(isnan(cal_SAL) & ~isnan(PSAL))='4';

  
    % ----- Add the flags from initial DMQC (in PREPARE_FLOATS) -----
    % Follow the rules in the manual and make the necessary variables.
    
    % 3.1. Editing raw qc flags in delayed-mode
    % If an error is identified, both PARAM_QC and PARAM_ADJUSTED_QC
    % should record ‘4’. Conversely, if good data have been identified
    % wrongly as bad by the real-time tests, then both PARAM_QC and
    % PARAM_ADJUSTED_QC should record ‘1’.
    % √ Done below in 'all parameters' section.
    
    % 3.2. Delayed-mode procedures for JULD, LATITUDE, LONGITUDE
    % JULD_QC = ‘8’.
    % POSITION_QC = ‘8’.
    % [] Not encountered and not implemented yet.
    
    % 3.3. Delayed-mode procedures for PRESSURE 
    % Bad data points identified by visual inspection from delayed-mode
    % analysts are recorded with PRES_ADJUSTED_QC = ‘4’ and PRES_QC =
    % '4'. For these bad data points, TEMP_QC, TEMP_ADJUSTED_QC, PSAL_QC,
    % PSAL_ADJUSTED_QC should also be set to ‘4’.
    % The direct flagging from RTQC of all parameters is done above.
    % E: In addition reversed PRESqcn give reversed T & S as well,
    % which is ok since their specific flags will be implemented below.
    PRES_QC=='4';                            TEMP_ADJUSTED_QC(ans)='4'; PSAL_ADJUSTED_QC(ans)='4'; % Old
    PRESqcn=='1'; PRES_ADJUSTED_QC(ans)='1'; TEMP_ADJUSTED_QC(ans)='1'; PSAL_ADJUSTED_QC(ans)='1'; % Reversed
    PRESqcn=='4'; PRES_ADJUSTED_QC(ans)='4'; TEMP_ADJUSTED_QC(ans)='4'; PSAL_ADJUSTED_QC(ans)='4'; % New
    % √ The point about PARAM_QC and FillValue is taken care of below directly in files.

    % Set pressure error and update scientific equation info accordingly:
    if ~isempty(PSAL_ADJUSTED_Cnew)
      PRES_ADJUSTED_ERROR = 2.5/6000 * PRES + 2;  
      'PRES_ADJUSTED_ERROR = 2.5/6000 * PRES + 2. ';
    else
      PRES_ADJUSTED_ERROR=ones(size(PRES_ADJUSTED))*2.4;  % also when SP is corrected for.
      'The quoted error is manufacturer specified accuracy.';
    end
    scientific_calib.comment.PRES(:,end+[1:size(ans,2)]) = repmat(ans,FN,1);

    % Please note that whenever PARAM_ADJUSTED_QC = ‘4’, both
    % PARAM_ADJUSTED and PARAM_ADJUSTED_ERROR should be set to FillValue.
    % √ Done below in the writing section.
    % [√] Also done for NaNs in PARAM_ADJUSTED.
  
    % 3.4. Delayed-mode procedures for TEMPERATURE 
    % Bad data points identified by visual inspection from delayed-mode
    % analysts are recorded with TEMP_ADJUSTED_QC = ‘4’ and TEMP_QC = '4'.
    TEMPqcn=='1'; TEMP_ADJUSTED_QC(ans)='1'; % Reversed
    TEMPqcn=='4'; TEMP_ADJUSTED_QC(ans)='4'; % New
    % √ The point about PARAM_QC is taken care of below.
  
    % √ TEMP_ADJUSTED, TEMP_ADJUSTED_ERROR, and TEMP_ADJUSTED_QC should be
    % filled even when the data are good and no adjustment is needed. In
    % these cases, TEMP_ADJUSTED_ERROR can be the manufacturer’s quoted
    % accuracy at deployment, which is 0.002°C.
    % Note: PARAM_ERROR does not exist!
    % √ Done on all three parameters below.
    TEMP_ADJUSTED_ERROR=ones(size(TEMP_ADJUSTED))*0.002; % Same for CP_cor adjustment.
  
    % Please use the SCIENTIFIC CALIBRATION section in the netCDF files to
    % record details of the delayed-mode adjustment.
    % [] 
  
    % 3.5. Delayed-mode procedures for SALINITY 
    % √ But we need to flag before OWC, so same procedure as TEMP is applied.
    PSALqcn=='1'; PSAL_ADJUSTED_QC(ans)='1'; % Reversed
    PSALqcn=='4'; PSAL_ADJUSTED_QC(ans)='4'; % New
    %if ~isempty(PSAL_ADJUSTED_Cnew)
    %if ~all(isnan((PSAL_ADJUSTED_Cnew)))
  
    if ~isempty(PSAL_ADJUSTED_Cnew)
      % PSAL_ADJUSTED_ERROR = minimum 0.004. DMQC operators should
      % consider increasing this error in case of sensor drift adjustment,
      % to take into account possible pressure dependency of the sensor
      % drift that cannot be detected with certainty.
      % [√] This is always set higher (at least 0.01) in the OWC
      % section above, in any case.
      % If no other error is found in delayed-mode and OWC salinity
      % adjustment is < 0.05, flag PSAL_QC='3' <= 2000 dbar:
      % [] ???
    end


  
    % ----- Changes to be done on all three parameters alike -----

    PARAM={'PRES','TEMP','PSAL'};  
    for j=1:length(PARAM)					% Loop parameters
  
      % Make sure the empty values do not have _ERROR (necessary after naive setting above) 
      eval([PARAM{j},'_ADJUSTED_ERROR(isnan(',PARAM{j},'_ADJUSTED))=NaN;'])
    
      % When all data in a profile is gone set qc='9' for adjusted 
      %eval([PARAM{j},'_QC(:,find(all(isnan(',PARAM{j},'))))=''9'';']) % set qc='9' for original 
      %eval([PARAM{j},'_ADJUSTED_QC(:,find(all(isnan(',PARAM{j},'_ADJUSTED))))=''9'';']) % set qc='4' for adjusted
      %eval([PARAM{j},'_ADJUSTED_QC(:,find(all(isnan(',PARAM{j},'_ADJUSTED) & ',PARAM{j},'_ADJUSTED_QC=='' '')))=''9'';']) %
      % BUT NOT IF THEY ARE 4 FOR BAD PROFILES SET IN PREVIOUS DMQC!
      % Old flags are already carried over, so if not changed back above, bad profiles should have all PSAL_QC==4 (and others).
      % So added to the test that all qc are not 4, before flagging with 9.
      % E.g., PSAL_ADJUSTED_QC(:, all NaNs & all not-fours & all empty new-qc) = '9'
      eval([PARAM{j},'_ADJUSTED_QC(:,find(all(isnan(',PARAM{j},'_ADJUSTED) & ',PARAM{j},'_QC~=''4'' & ',PARAM{j},'_ADJUSTED_QC=='' '')))=''9'';']) %
    
      % Make sure flags (currently only '4' and '1', and '9' added) match
      % between the two 'levels':
      eval([PARAM{j},'_QC=',PARAM{j},'_ADJUSTED_QC;']) 

      % LAST but not least:
      % The variable PROFILE_<PARAM>_QC should be recomputed when <PARAM>_ADJUSTED_QC becomes available.
      % N is defined as the percentage of levels with good data where:
      % • QC flag values of 1, 2, 5, or 8 are GOOD data
      % • QC flag values of 9 (missing) are NOT USED in the computation
      % • All other QC flag values are BAD data
      % The computation should be taken from <PARAM>_ADJUSTED_QC if available and from
      % <PARAM>_QC otherwise.
      % n Meaning
      % “ “ No QC was performed
      % A N = 100%; All profile levels contain good data
      % B 75% <= N < 100%
      % C 50% <= N < 75%
      % D 25% <= N < 50%
      % E 0% < N < 25%
      % F N = 0%; No profile levels have good data
      eval(['PROFILE_',PARAM{j},'_QC=repmat(char(32),1,size(',PARAM{j},'_ADJUSTED_QC,2));']);
      eval(['sum(',PARAM{j},'_ADJUSTED_QC==''1'' | ',PARAM{j},'_ADJUSTED_QC==''2'' | ',PARAM{j},'_ADJUSTED_QC==''5'' | ',PARAM{j},'_ADJUSTED_QC==''8'')./PROFILE_',PARAM{j},'_N*100;']);
      eval(['PROFILE_',PARAM{j},'_QC(ans==0) =''F''; PROFILE_',PARAM{j},'_QC(ans>0)  =''E''; PROFILE_',PARAM{j},'_QC(ans>=25) =''D'';']);
      eval(['PROFILE_',PARAM{j},'_QC(ans>=50)=''C''; PROFILE_',PARAM{j},'_QC(ans>=75)=''B''; PROFILE_',PARAM{j},'_QC(ans==100)=''A'';']);
  
    end

    % ----- Flag the deep corrected salinities -----
    % 3.10.1 Summary of delayed-mode QC flag scheme for Deep-Argo data
    % Since PSAL_ADJUSTED_QC and PSAL_QC will be different in ths case,
    % this must be done after the common section above.
    if ~isempty(PSAL_ADJUSTED_Cnew)
      % If no other error is found in delayed-mode and OWC salinity
      % adjustment is < 0.05, flag PSAL_QC='3' > 2000 dbar:
      PRES_ADJUSTED>2000 & PSAL_ADJUSTED_QC=='1'; PSAL_QC(ans)='3';
    end


    % ------ FINALLY WRITE TO THE D-FILES -----
  
    for i=1:FN							% Loop files, i.e., columns
    
      ncid = netcdf.open(Dfiles{i},'WRITE');			% Open D-file for writing
      dimid = netcdf.inqDimID(ncid,'N_LEVELS'); [~,M] = netcdf.inqDim(ncid,dimid);% Get length of profile in file
      dimid = netcdf.inqDimID(ncid,'N_PROF'); [~,N] = netcdf.inqDim(ncid,dimid);	% Get number of profiles in file

      for j=1:length(PARAM)					% Loop parameters
  
	% Update PARAM_ADJUSTED_QC:
	adjqcid = netcdf.inqVarID(ncid,[PARAM{j},'_ADJUSTED_QC']);% Get adjusted parameter qc ID
	adjqc=netcdf.getVar(ncid,adjqcid);			% Read adjusted parameter qc
	eval(['adjqc(:,1)=',PARAM{j},'_ADJUSTED_QC(1:M,i);']);	% Overwrite with PARAM_ADJUSTED_QC
	netcdf.putVar(ncid,adjqcid,adjqc);			% Write back
	
	% Update PARAM_QC: 
	qcid = netcdf.inqVarID(ncid,[PARAM{j},'_QC']);		% Get parameter qc ID
	qc=netcdf.getVar(ncid,qcid);				% Read parameter qc
	eval(['qc(:,1)=',PARAM{j},'_QC(1:M,i);']);		% Overwrite with PARAM_QC
	netcdf.putVar(ncid,qcid,qc);				% Write back

	% Update PROFILE_PARAM_QC: 
	pqcid = netcdf.inqVarID(ncid,['PROFILE_',PARAM{j},'_QC']);% Get profile parameter qc ID
	pqc=netcdf.getVar(ncid,pqcid);				% Read profile parameter qc
	eval(['pqc(1,1)=PROFILE_',PARAM{j},'_QC(1,i);']);		% Overwrite with PROFILE_PARAM_QC
	netcdf.putVar(ncid,pqcid,pqc);				% Write back
      
	% Update PARAM_ADJUSTED: 
	adjid = netcdf.inqVarID(ncid,[PARAM{j},'_ADJUSTED']);	% Get adjusted parameter ID
	adj=netcdf.getVar(ncid,adjid);				% Read adjusted parameter
	eval(['adj(:,1)=',PARAM{j},'_ADJUSTED(1:M,i);']);		% Overwrite with PARAM_ADJUSTED

	% Update PARAM_ADJUSTED_ERROR:
	errid = netcdf.inqVarID(ncid,[PARAM{j},'_ADJUSTED_ERROR']); % Get adjusted parameter err ID
	err=netcdf.getVar(ncid,errid);				% Read adjusted parameter err
	eval(['err(:,1)=',PARAM{j},'_ADJUSTED_ERROR(1:M,i);']);	% Overwrite with PARAM_ADJUSTED_ERROR

	% Prepare for porting PARAM_* to PARAM_ADJUSTED_*:
	% Fill with fillValue only if needed:
	%eval([PARAM{j},'_ADJUSTED_QC(1:M,i)==''4'' | isnan(',PARAM{j},'_ADJUSTED_QC(1:M,i));']);	% Check for new flags ('4'):
	eval([PARAM{j},'_ADJUSTED_QC(1:M,i)==''4'' | isnan(',PARAM{j},'_ADJUSTED(1:M,i));']);	% Check for new flags ('4') or NaNs in data
	if any(ans,'all')						% If there are any,
	  % Put fillvalues in PARAM_ADJUSTED:  
	  [~,fillValue] = netcdf.inqVarFill(ncid,adjid);		% get PARAM_ADJUSTED's fillvalue
	  adj(find(ans),1)=fillValue;				% fill the data from PARAM. 
          % Put fillvalues in PARAM_ADJUSTED_ERROR: 
	  [~,fillValue] = netcdf.inqVarFill(ncid,errid);		% get PARAM_ADJUSTED_ERROR's fillvalue
	  err(find(ans),1)=fillValue;				% fill err from PARAM. 
	end
	% netcdf.putVar(ncid,adjid,adj);				% Write PARAM over to PARAM_ADJUSTED (in any case)
	% netcdf.putVar(ncid,errid,err);				% Write PARAM_ERROR over to PARAM_ADJUSTED_ERROR (in any case)
	netcdf.putVar(ncid,adjid,adj);				% Write back to PARAM_ADJUSTED
	netcdf.putVar(ncid,errid,err);				% Write back to PARAM_ADJUSTED_ERROR
    
      end	% Loop parameters

      % The variable DATA_MODE should record ‘D’.
      dmid = netcdf.inqVarID(ncid,'DATA_MODE');			% Get ID
      dm=netcdf.getVar(ncid,dmid);				% Read 
      dm(1)='D';							% Overwrite 
      netcdf.putVar(ncid,dmid,dm);				% Write back 

      % The variable DATA_STATE_INDICATOR should record ‘2C’ or ‘2C+’.
      dsid = netcdf.inqVarID(ncid,'DATA_STATE_INDICATOR');	% Get ID
      ds=netcdf.getVar(ncid,dsid);				% Read 
      '2C'; ds(1:2,1)=ans';					% Overwrite 
      netcdf.putVar(ncid,dsid,ds);				% Write back 

      % The variable DATE_UPDATE should record the date of last update of the netcdf file, in the format YYYYMMDDHHMISS.
      daid = netcdf.inqVarID(ncid,'DATE_UPDATE');			% Get ID
      da=replace(datestr(now,30),'T','');				% Make date (also used below)
      netcdf.putVar(ncid,daid,da);				% Write date

      % GLOBAL ATTRIBUTES
      % Close file and start using the simple to use functions, because
      % of some trouble with the low-level functions.
      netcdf.close(ncid); 
    
      % A history record should be appended to the HISTORY section of the netcdf file:
      hi=ncreadatt(Dfiles{i},'/','history');				% Read
      hi=[hi,'; ',datestr2(now,111),'Z ',history_update]; % Add to string
      ncwriteatt(Dfiles{i},'/','history',hi);				% Write back 
    
      % Update the user-manual version:
      ncwriteatt(Dfiles{i},'/','user_manual_version',user_manual_version);		% Write 
    
      % DMQC operator:
      %comment_dmqc_operator1 = "PRIMARY | http://orcid.org/16-digit-number 1 | operator name 1, institution 1"
      if length(dmqc_operator)==1
	'comment_dmqc_operator';
	ncwriteatt(Dfiles{i},'/',ans,[dmqc_operator.parameter,' | http://orcid.org/',dmqc_operator.orcid,' | ',dmqc_operator.name,', ',dmqc_operator.institution]);
      else
	for j=1:length(dmqc_operator)
	  ['comment_dmqc_operator',int2str(j)];
	  ncwriteatt(Dfiles{i},'/',ans,[dmqc_operator(j).parameter,' | http://orcid.org/',dmqc_operator(j).orcid,' | ',dmqc_operator(j).name,', ',dmqc_operator(j).institution]);
	end
      end

    
      % 3.6.2. Scientific calibration information for each profile
      % Read what's in file already:
      ncread(Dfiles{i},'PARAMETER'); vars=ans(:,:,1,1)';
      % Check if all vars in this particular file have their scientific_calib info:
      if ~all(ismember(sort(strip(string(vars))),sort(string(fieldnames(scientific_calib.equation)))))
	error(['List in ''PARAMETER'' in ',Dfiles{i},...
	       ' does not have all scientific_calib struct from tests!']); 
      end
      % if size(vars,1)>3, 
      %   warning(['More than 3 variables in ',Dfiles{i}]); 
      %   vars=vars(contains(string(vars),{'PRES','TEMP','PSAL'}),:);
      % end
      calinfo={'EQUATION','COEFFICIENT','COMMENT','DATE'};
      for k=1:length(calinfo) % Loop calib struct fields
	scal=ncread(Rfiles{i},['SCIENTIFIC_CALIB_',calinfo{k}]); 
	D=size(scal);
	SCAL=repmat(char(32),D); % New matrix 
	if length(D)>=4 & D(4)>1, SCAL(:,:,:,2:end)=scal(:,:,:,2:end); end
	scal=scal(:,:,1,1)'; % Change only the first profile of the two
	for j=1:size(vars,1) % Loop variables
	  % Combine old and new information:
	  scal(j,:); % from R-file
	  replace(ans,'none','');replace(ans,char(0),char(32));strip(ans);if isempty(ans),ans(1:4)='none';end % Remove superfluous 'none'
	  oldinfo=ans;
	  scientific_calib.(lower(calinfo{k})).(strip(vars(j,:)))(i,:);  % from your current DMQC (e.g., scientific_calib.equation.PRES(i,:))
	  replace(ans,'none','');replace(ans,char(0),char(32));strip(ans);if isempty(ans),ans(1:4)='none';end % Remove superfluous 'none'
	  newinfo=ans;
	  if strcmp(calinfo{k},'DATE') | strcmp(oldinfo,newinfo) | contains(oldinfo,'none')  % just overwrite
	    ans=[newinfo];
	  else
	    ans=[newinfo];
	    %ans=[oldinfo,'. ',newinfo]; % No, we do not carry old scientific info from RTQC through DMQC!
	  end
	  % Clean up the new string:
	  replace(ans,'none','');replace(ans,char(0),char(32));strip(ans);if isempty(ans),ans(1:4)='none';end % Remove superfluous 'none' and change char(0)
	  replace(ans,{'='},{'= '}); replace(ans,{'+/-'},{'+/- '});			% Fix the issue of disappearing spaces
	  replace(ans,{'  '},{' '});							% But remove any double spaces
	  replace(ans,{'    '},{''});							% After that remove large spaces completely (big gaps)
	  scal(j,:)=' ';									% Clean line in scal for the file	
	  scal(j,1:length(ans))=ans;							% Update line in scal for the file	
	  scientific_calib.(lower(calinfo{k})).(strip(vars(j,:)))(i,1:end)=' ';		% Clean line in the internal object 
	  scientific_calib.(lower(calinfo{k})).(strip(vars(j,:)))(i,1:size(ans,2))=ans;	% Update line in the internal object 
	  if length(ans)>256, error(['SCIENTIFIC_CALIB_',calinfo{k},'_',strip(vars(j,:)),' = ''',ans,''' is longer than 256!']); end
	end
	SCAL(:,:,1,1)=scal';								% Only on the first N_PROF
	ncwrite(Dfiles{i},['SCIENTIFIC_CALIB_',calinfo{k}],SCAL);
	%ncread(Dfiles{i},['SCIENTIFIC_CALIB_',calinfo{k}]); ans(:,:,1,1)'
      end
      % There should be stucture for each, similar to
      % scientific_calib.equation.PRES;
      % scientific_calib.equation.TEMP;
      % scientific_calib.equation.PSAL;
      % for the parameters that exists in file. If file has more
      % parameters than made in PREPARE_FLOATS, an error will occur,
      % and you will have to update PREPARE_FLOATS to include those.
      %
      % [] These might have to be made for each cycle, later.
    

    end	% Loop files
  
  
    % ------- Write the scientific calibration to a tex file for report appendix --------
    disp(['Writing the scientific calibration info to report appendix for float ',float_names{I},'.']);
    vars=fieldnames(scientific_calib.comment);
    fid=fopen('scientific_calib_tabular.tex','w');
    fprintf(fid,'%s\n',['\subsection*{Scientific calibration information}']);
    fprintf(fid,'%s\n',['The scientific calibration information written to the D-files are summarized in',...
			' Tables~\ref{tab:scientific_calib_',vars{1},'}--\ref{tab:scientific_calib_',vars{end},'}.']);
    fprintf(fid,'%s\n',['Note that adjustements are done and registered regardless of the amount of valid', ...
			' data, hence some cycles will record scientific calibration information even when the whole profile consists of fillvalues.']);
    for k=1:length(vars) % Loop vars
      linecount=1;
      fprintf(fid,'%s\n','\begin{table}[!h]');
      fprintf(fid,'%s%s%s%s%s%s\n','\caption{Information filled in the SCIENTIFIC\_CALIB section for ',replace(vars{k},'_','\_'),' in the D-files. Empty string means parameter not present in file.}',...
	      ' \label{tab:scientific_calib_',vars{k},'}');
      %fprintf(fid,'%s\n','\caption{Information filled in the SCIENTIFIC\_CALIB section for the variables, in the D-files.} \label{tab:scientific_calib}');
      fprintf(fid,'%s\n','\begin{tabular}[b]{|r|l|p{3cm}|p{9cm}|}');
      fprintf(fid,'%s\n','\hline');
      fprintf(fid,'%s & %s & %s & %s\\\\ \n','Parameter','Field','Cycles/files','String');
      fprintf(fid,'%s\n','\hline');
      fprintf(fid,'%s\n','\hline');
      fprintf(fid,'%s & ',replace(vars{k},'_','\_'));
      for j=1:length(calinfo) % Loop calib struct fields
	if j==1, form='%s & '; else, form='& %s & '; end
	fprintf(fid,form,calinfo{j});
	%if j==2 & k==3, return; end
	[C,IA,IC]=unique(cellstr(scientific_calib.(lower(calinfo{j})).(vars{k})),'stable');
	for i=1:length(C)
	  if i==1, form='%s & %s\\\\ \n'; else, form='& & %s & %s\\\\ \n'; end
	  replace(C{i},'_','\_');						% Replace for LaTeX
	  replace(ans,{'='},{'= '}); replace(ans,{'+/-'},{'+/- '});	% Fix the issue of disappearing spaces
	  replace(ans,{'  '},{' '});					% remove any double spaces
	  lans=length(ans);
	  %if length(find(IC==i))==2, return; end
	  fprintf(fid,form,zipnumstr(CYCLE_NUMBER(find(IC==i)),', '),ans);
	  %fprintf(fid,form,zipnumstr(find(IC==i),', '),ans);
	  linecount=linecount+ceil(lans/40);	% Increase by approximate line number
	  if linecount>52				% Split table at reasonable length
	    % end current table
	    fprintf(fid,'%s\n','\hline');
	    fprintf(fid,'%s\n','\end{tabular}');
	    fprintf(fid,'%s\n','\end{table}');
	    % begin new table 
	    fprintf(fid,'%s\n','\addtocounter{table}{-1}');
	    fprintf(fid,'%s\n','\begin{table}[!p]');
	    fprintf(fid,'%s\n','\caption{Continued.}');
	    %fprintf(fid,'%s\n','\caption{Table~\protect\ref{tab:scientific_calib} continued.} \label{}');
	    fprintf(fid,'%s\n','\begin{tabular}[b]{|r|l|p{3cm}|p{9cm}|}');
	    fprintf(fid,'%s\n','\hline');
	    fprintf(fid,'%s & %s & %s & %s\\\\ \n','Parameter','Field','Cycles/files','String');
	    fprintf(fid,'%s\n','\hline');
	    fprintf(fid,'%s\n','\hline');
	    linecount=1;
	  end
	end
      end % Loop calib struct fields
      fprintf(fid,'%s\n','\hline');
      fprintf(fid,'%s\n','\end{tabular}');
      fprintf(fid,'%s\n','\end{table}');
    end % Loop vars 
    fclose(fid);

  
  
  end	% Loop floats

end	% Loop direction


cd(my_working_dir); % Go back to the main working dir

% Useful lines that can be inserted above for debugging purposes:
% strcat(PSALqcn(1:2,:),'|',PSAL_QC(1:2,:)) % early test 
% strcat(PSALqcn(1:2,:),'|',PSAL_QC(1:2,:),'|',PSAL_ADJUSTED_QC(1:2,:)) % test



