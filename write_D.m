% This script writes the results of DMQC to D-files.

clear all; close all
init_dmqc; % Paths and filenames.

for I=1:length(download_dir)	% Loop floats

  load(outfiles{I}); 
      % 'LAT','LONG','DATES','PRES','SAL','TEMP','PTMP','PROFILE_NO',...
      % 'JULD','CYCLE_NUMBER',...
      % 'PROFILE_*_N','*qco','*qcn','Rfiles','scientific_calib','PSAL_ADJUSTED_Cnew',...
      % 'ctdmodel','PSAL'); 

  %PSAL=SAL; 
  clear SAL
  load([calib_dir,'cal_',float_names{I}]); % results of OWC
  load([calib_dir,'calseries_',float_names{I}]); % results of OWC
  cd([my_working_dir,'DMQC',filesep,float_names{I}]); % Just to be safe

  % Copy and rename Rfiles to Dfiles:
  if exist(rootdirout{I},'dir'), rmdir(rootdirout{I},'s'), end
  mkdir(rootdirout{I});  
  Dfiles=replace(replace(Rfiles,'/R','/D'),rootdirin{I},rootdirout{I});
  command=strcat({'cp '},Rfiles,{' '},Dfiles);
  for i=1:length(command), system(command{i}); end

  % Count and test files vs. matrix data:
  FN=length(Dfiles);
  if FN~=size(PRES,2), error('Number of files does not match number of columns in data!'); end
  
  % Initiate variables for NetCDF (PSAL below):
  PRES_ADJUSTED=PRES;		TEMP_ADJUSTED=TEMP;  
  PRES_ADJUSTED_QC=PRESqco;	TEMP_ADJUSTED_QC=TEMPqco; % Carry over the old flags
  %clear PSAL_ADJUSTED*
  PRES_QC=PRESqco; TEMP_QC=TEMPqco; PSAL_QC=SALqco; % Rename old flags to fit netcdf names

  % Update scientific equation info accordingly:
  'PRES_ADJUSTED = PRES. '; %'none';
	scientific_calib.equation.PRES(:,end+[1:size(ans,2)]) = repmat(ans,FN,1);
  'TEMP_ADJUSTED = TEMP. '; %'none';
	scientific_calib.equation.TEMP(:,end+[1:size(ans,2)]) = repmat(ans,FN,1);

  
  % SAL - from the pressure dep bias. if strcmp(platyp,'ARVOR-D');
  % must be put as PSAL_ADJUSTED_Cnew. Use PSAL, since it has the
  % flagged data NaNed out.
  
  
  % ----- ADD THE RESULTS OF OWC TO OBJECTS -----------------------
  [PSAL_ADJUSTED,PSAL_ADJUSTED_ERROR]=deal(nan(size(PSAL))); 
  PSAL_ADJUSTED_QC=repmat(char(32),size(PSAL));
  Na=length(cal_action{I});
  if length(unique(calseries(calseries>0)))~=Na, error('Number of specified actions differ from parts of calseries!'); end
  for i=1:Na			% Loop actions
    jj=find(calseries==i);	% True at profiles for i-th action 
    switch cal_action{I}(i);	% What is i-th action?
     case 0			% Data are good; no adjustment has been applied:
      PSAL_ADJUSTED_QC(:,jj) = PSAL_QC(:,jj); % RTQC flags (DMQC flags are added below)
      PSAL_ADJUSTED_ERROR(:,jj) = max(cal_SAL_err(:,jj),0.01); % Yes, also when no adjustment is done.
      if isempty(PSAL_ADJUSTED_Cnew)
	PSAL_ADJUSTED(:,jj) = PSAL(:,jj); % Original PSAL with NaNs assigned by PREPARE_FLOATS
	'PSAL_ADJUSTED = PSAL. '; %'none';
		scientific_calib.equation.PSAL(jj,end+[1:size(ans,2)]) = repmat(ans,length(jj),1);
      else
	PSAL_ADJUSTED(:,jj) = PSAL_ADJUSTED_Cnew(:,jj); % Original PSAL with NaNs assigned by PREPARE_FLOATS
	'PSAL_ADJUSTED = PSAL_ADJUSTED_Cnew. '; %'none';
		scientific_calib.equation.PSAL(jj,end+[1:size(ans,2)]) = repmat(ans,length(jj),1);
      end
      'none';	
		scientific_calib.coefficient.PSAL(jj,end+[1:size(ans,2)]) = repmat(ans,length(jj),1);
      'No significant salinity offset or drift detected. The quoted error is max[0.01, statistical uncertainty] in PSS-78. ';
		scientific_calib.comment.PSAL(jj,end+[1:size(ans,2)]) = repmat(ans,length(jj),1);
      %W PSAL_ADJUSTED = PSAL (original value);
      %W PSAL_ADJUSTED_QC = ‘1’ or ‘2’;
      %W PSAL_ADJUSTED_ERROR = maximum [statistical uncertainty, 0.01].  
     case 1			% Data show sensor drift or offset; adjustment has been applied:
      PSAL_ADJUSTED(:,jj) = cal_SAL(:,jj);     % OWC carries over the CPnew corrections together with its adjustment, in cal_SAL
					       % But the scientifc cal equation becomes different.
      PSAL_ADJUSTED_QC(:,jj) = PSAL_QC(:,jj); % RTQC (DMQC flags are added below)
      PSAL_ADJUSTED_ERROR(:,jj) = max(cal_SAL_err(:,jj),0.01);
      if isempty(PSAL_ADJUSTED_Cnew)
	'PSAL_ADJUSTED = PSAL + dS. ';
		scientific_calib.equation.PSAL(jj,end+[1:size(ans,2)]) = repmat(ans,length(jj),1);
      else
	'PSAL_ADJUSTED = PSAL_ADJUSTED_Cnew + dS. ';
		scientific_calib.equation.PSAL(jj,end+[1:size(ans,2)]) = repmat(ans,length(jj),1);
      end
      % 'Vertically averaged dS = 0.012 +/- 0.003';
      %strcat('Vertically averaged dS = ',num2str(median([cal_SAL(:,jj)-PSAL(:,jj)],1,'omitnan')','%5.3f'),' +/- ',num2str(median(cal_SAL_err(:,jj),1,'omitnan')','%5.3f'));
      strcat('Vertically averaged dS = ',num2str(mean([cal_SAL(:,jj)-PSAL(:,jj)],1,'omitnan')','%5.3f'),' +/- ',num2str(mean(cal_SAL_err(:,jj),1,'omitnan')','%5.3f'),'. ');
		scientific_calib.coefficient.PSAL(jj,end+[1:size(ans,2)]) = ans;
      'Significant salinity sensor offset or drift detected. OW least squares fit adopted. The quoted error is max[0.01, statistical uncertainty] in PSS-78. ';
		scientific_calib.comment.PSAL(jj,end+[1:size(ans,2)]) = repmat(ans,length(jj),1);
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
      %W "PSAL_ADJUSTED = FillValue;"			
      %W "PSAL_ADJUSTED_QC = ‘4’; PSAL_QC = '4';"	
      %W "PSAL_ADJUSTED_ERROR = FillValue."
      %W   scientific_calib.equation.PSAL = 'none';
      %W   scientific_calib.coefficient.PSAL = 'none';
      %W   scientific_calib.comment.PSAL = 'Salinity data are bad and unadjustable because of the Druck oil microleak problem.';
      % NaNs are replaced by FillValues below.
    end % What i-th action switch
  end % Loop actions
  % This naive adding to the scientific calib info, often results in
  % superfluous 'none' in the texts. This will be cleaned up below (3.6.2).
  scientific_calib.date.PSAL  	= repmat(replace(datestr(now,30),'T',''),FN,1); % Always replace

  % Some times OWC removes data above a set of NaN values; flag these
  % as bad (qc='4'):
  SALqcn(isnan(cal_SAL) & ~isnan(PSAL))=4;
  
  
  % ------ ADD THE FLAGS FROM INITIAL DMQC (PREPARE_FLOATS) ---------------
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
  PRESqcn==1; PRES_ADJUSTED_QC(ans)='1'; TEMP_ADJUSTED_QC(ans)='1'; PSAL_ADJUSTED_QC(ans)='1';
  PRESqcn==4; PRES_ADJUSTED_QC(ans)='4'; TEMP_ADJUSTED_QC(ans)='4'; PSAL_ADJUSTED_QC(ans)='4';
  PRESqco=='4';                          TEMP_ADJUSTED_QC(ans)='4'; PSAL_ADJUSTED_QC(ans)='4';
  % √ The point about PARAM_QC and FillValue is taken care of below directly in files.
  clear *qco

  % Set pressure error and update scientific equation info accordingly:
  if ~isempty(PSAL_ADJUSTED_Cnew)
    PRES_ADJUSTED_ERROR = 2.5/6000 * PRES + 2;  
    'PRES_ADJUSTED_ERROR = 2.5/6000 * PRES + 2. ';
  else
    PRES_ADJUSTED_ERROR=ones(size(PRES_ADJUSTED))*2.4; 
    'PRES_ADJUSTED_ERROR = 2.4. ';
  end
  scientific_calib.comment.PRES(:,end+[1:size(ans,2)]) = repmat(ans,FN,1);

  % Please note that whenever PARAM_ADJUSTED_QC = ‘4’, both
  % PARAM_ADJUSTED and PARAM_ADJUSTED_ERROR should be set to FillValue.
  % √ Done below in the writing section.
  % [√] Also done for NaNs in PARAM_ADJUSTED.
  
  % 3.4. Delayed-mode procedures for TEMPERATURE 
  % Bad data points identified by visual inspection from delayed-mode
  % analysts are recorded with TEMP_ADJUSTED_QC = ‘4’ and TEMP_QC = '4'.
  TEMPqcn==1; TEMP_ADJUSTED_QC(ans)='1';
  TEMPqcn==4; TEMP_ADJUSTED_QC(ans)='4';
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
  SALqcn==1; PSAL_ADJUSTED_QC(ans)='1';
  SALqcn==4; PSAL_ADJUSTED_QC(ans)='4';
  % if ~isempty(PSAL_ADJUSTED_Cnew)
  %   % PSAL_ADJUSTED_ERROR = minimum 0.004. DMQC operators should
  %   % consider increasing this error in case of sensor drift adjustment,
  %   % to take into account possible pressure dependency of the sensor
  %   % drift that cannot be detected with certainty.
  %   % [√] This is always set higher (at least 0.01) in the OWC
  %   % section above, in any case.
  % end

  % ------ CHANGES TO BE DONE ON ALL THREE PARAMETERS ALIKE -------------------

  PARAM={'PRES','TEMP','PSAL'};  
  for j=1:length(PARAM)					% Loop parameters
  
    % Make sure the empty values do not have _ERROR (necessary after naive setting above) 
    eval([PARAM{j},'_ADJUSTED_ERROR(isnan(',PARAM{j},'_ADJUSTED))=NaN;'])
    
    % When all data in a profile is gone set qc='9' for adjusted
    %eval([PARAM{j},'_QC(:,find(all(isnan(',PARAM{j},'))))=''9'';']) % set qc='9' for original 
    %eval([PARAM{j},'_ADJUSTED_QC(:,find(all(isnan(',PARAM{j},'_ADJUSTED))))=''9'';']) % set qc='4' for adjusted
    eval([PARAM{j},'_ADJUSTED_QC(:,find(all(isnan(',PARAM{j},'_ADJUSTED) & ',PARAM{j},'_ADJUSTED_QC=='' '')))=''9'';']) %
     
    % make sure flags (currently only '4' and '1', and '9' added) match
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
  
  
  
  % ------ WRITE TO THE FILES --------------------------------------------
  
  for i=1:FN							% Loop files, i.e., columns
    ncid = netcdf.open(Dfiles{i},'WRITE');
    dimid = netcdf.inqDimID(ncid,'N_LEVELS'); [~,M] = netcdf.inqDim(ncid,dimid);% Get length of profile in file
    dimid = netcdf.inqDimID(ncid,'N_PROF'); [~,N] = netcdf.inqDim(ncid,dimid);	% Get number of profiles in file
    PARAM={'PRES','TEMP','PSAL'};  
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
      % varid = netcdf.inqVarID(ncid,PARAM{j});			% Get parameter ID
      % data=netcdf.getVar(ncid,varid);				% Read parameter
      % adjid = netcdf.inqVarID(ncid,[PARAM{j},'_ADJUSTED']);	% Get adjusted parameter ID
      % errid = netcdf.inqVarID(ncid,[PARAM{j},'_ADJUSTED_ERROR']);% Get parameter err ID
      % err=netcdf.getVar(ncid,errid);				% Read parameter err
      % if strcmp(PARAM{j},'TEMP'), err(:,:)=0.002; end		% Fill the err with set value
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
    if size(vars,1)>3, warning(['More than 3 variables in ',Dfiles{i}]); end
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
  fid=fopen('scientific_calib_tabular.tex','w');
  linecount=1;
  fprintf(fid,'%s\n','\begin{table}[!h]');
  fprintf(fid,'%s\n','\caption{Information filled in the SCIENTIFIC\_CALIB section for the variables, in the D-files.} \label{tab:scientific_calib}');
  fprintf(fid,'%s\n','\begin{tabular}[b]{|r|l|p{3cm}|p{9cm}|}');
  fprintf(fid,'%s\n','\hline');
  fprintf(fid,'%s & %s & %s & %s\\\\ \n','Parameter','Field','Cycles/files','String');
  fprintf(fid,'%s\n','\hline');
  fprintf(fid,'%s\n','\hline');
  vars=fieldnames(scientific_calib.comment);
  for k=1:length(vars)
    %fprintf(fid,'%s & ',['Scientific calibration information for ',vars{k},':']);
    fprintf(fid,'%s & ',vars{k});
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
	fprintf(fid,form,zipnumstr(find(IC==i),', '),ans);
	linecount=linecount+ceil(lans/40);	% Increase by approximate line number
	if linecount>60				% Split table at reasonable length
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
    end
    fprintf(fid,'%s\n','\hline');
  end
  %fprintf(fid,'%s\n','\hline');
  fprintf(fid,'%s\n','\end{tabular}');
  fprintf(fid,'%s\n','\end{table}');
  fclose(fid);

end	% Loop floats


