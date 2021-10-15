% INIT_DMQC.M is init-file and setup instructions.  
% DMQC-fun v0.9.1; See Contents.m for overview; work_log.txt for workflow description.
%%%%%%%%%% Read further down below for instructions on how to set up your system! %%%%%%%%%%%%%%%%  
%%%%%%%%%% The first part is where you control everything! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------ LIST ALL YOUR FLOATS HERE: ----------------------------
TCL=complex(0,1); % Do not touch!

% index =     1         2         3         4         5         6         7         8         9         10        11        12        13        14        15        16        17        18*       19*       20        21        22        23        24        25
float_names={'6903549','6903550','6903551','6903574','6903552','6903553','6903554','6903555','6903567','6903568','6903569','6903570','6903559','6903560','6903561','6903562','6903563','6903564','6903565','6903566','6903556','6903557','6903558','6903571','6903573'};
cal_action ={ [0]     , [0]     , [1]     , [1]     , [0]     , [0]     , [0]     , [0]     , [0]     , [0]     , [0]     , [0]     , [0]     , [0]     , [0 1]   , [0 1 4] , [1]     , [0]     , [0]     , [0]     ,  [0]    , [0]     , [0]     , [0]     , [0]      }; 
checked    ={ 177       147       124       80        41        113       152       70        60        89        115       82        85        86        86        85        86        156       101       51      , 30        67        66        31        20      };
Nclu       ={ 2       2           2         2         2         2         2         2         2         2         2         2         2         2         3*TCL     2*TCL     1*TCL     2         2         2         2         2         2         2         2       };

% ------------ FLOAT SELECTION: ----------------------------------------
% You can select which floats to analyse by just giving specific index
% numbers and all the parameters will be reduced to your selection. For
% example like this, and uncomment the subset to DMQC:
% 1:length(float_names) % All floats
1		% Just the first, as default.
% 1:4 		% BGC (PROVOR CTS4)
% [5:8 10:12]	% Bio (APEX AFP11) 
% [13:17 20]  	% Core (Arvor)
% 21:25  	% Dyp (Arvor Deep) 
% 18:19  	% Those two floats in Barents Sea, no ref data, no OWC.
% 5		% The first APEX float
% 12		% The last APEX, on greylist
% [1:17 20:25] % All but those two that cannot be calibrated due to lack of refdata. 
% 23		% A deep arvor to test on
float_names=float_names(ans); Nclu=Nclu(ans); cal_action=cal_action(ans); checked=checked(ans); 

% index is just a reference for your eyeing of the 'table' of floats
%	and parameters, and thus that line can stay commented. 
% float_names is necessary for ow_calibration.m. DO NOT change the
%	name of this variable!
% cal_action is your calibration decision. Set this initially to [0]
%	for every float.
% checked is the mumber of profiles the RTQC flags is already checked for
%	(you will be prompted to update this number in the
%	process). Set this initially to 0 for every float. 
% Nclu is number of clusters to use in the visual comparison with
%	reference data. Set to 2 or longer if float has traversed
%	several basins or regions. Multiply with TCL if time (i.e.,
%	cycle number) is to affect the clustering. 


% ------------- OVERVIEW OF THIS INIT FILE: --------------------------
%
% This file will be used by all the main DMQC-fun scripts, as well as
% OWC functions, when using the DMQC-fun setup. It has the following
% sections:
% - A list of all your floats and corresponding parameters you will
%   change during DMQC (above for convenience). 
% - Float selection part for which floats to DMQC (above).
% - This overview of this file.
% - First-time preparations, i.e., the install instructions.
% - Operator name etc. for the report and D-files.
% - Assorted parameters.
% - Paths for everything:
%	MATLABOW directories, according to your installation of OWC toolbox;
%	Download directories for float data (your choice);
%	Reference data directories (your choice);
%	List of WMO-squares for yor area of interest (your choice);
%	Your other relevant Argo directories (your choice);
%	Putting paths in order (fixed);
%	Auxillary TOOLBOXES included with DMQC-fun (fixed);
%	File names for PREPARE_PROFILES (fixed);
%	Directories for WRITE_D (fixed);
% - Finally, the first build of your working directories, a part of
%   the script that copies the necessary files into separate folders
%   for each of your floats (see next section).


% ------------- FIRST-TIME SETUP AND PREPARATIONS: -------------------

% 1) Install the following toolboxes in a designated directory for your
% Argo related toolboxes (e.g., /Users/a21627/matlab/toolbox/Argo/):
% 
%	DMQC-fun - https://github.com/imab4bsh/DMQC-fun.git
%	matlab_owc -  https://github.com/ArgoDMQC/matlabow.git
% 
% Also install 
%
%	evenmat - https://github.com/evenrev1/evenmat.git
%	m_map - https://www.eoas.ubc.ca/~rich/map.html
%	GSW - http://www.teos-10.org/pubs/gsw/html/
%
% The oceans, seawater, and some other necessary toolboxes are provided
% herein. Unzip and add their paths if you don't have them installed
% already.
%
% When adding the paths to the toolboxes, good practice is to make
% sure DMQC-fun comes above these other toolboxes in the path list
% (i.e., addpath it last), so that our altered versions of some files
% (e.g., plot_diagnostics_ow) are used instead of the matlab_owc
% version. However, there is a fail-safe mechanism buildt into INIT_DMQC
% that adds and orders these paths upon first time use.
% 
% LaTeX: DMQC-fun includes a LaTeX report template and the Matlab
% scripts produce snippets of content linked into that template. You
% will need a working version of LaTeX, and for instance dvipdfm to make
% PDF. However, this is not necessary for DMQC-fun to be useful. (You
% may even be able to link the produced figures and text parts into some
% other word processor of choice.)

% 2) You need to copy some files from the toolboxes manually, as you do
% not want updating of toolboxes to erase your own settings. From
% DMQC-fun toolbox copy the following files to your working directory
% for Argo DMQC (my_working_dir):
%
%	init_dmqc.m (this file) 
%	work_log.txt
%
% A path to your working directory will be added to make sure this
% local INIT_DMQC is the one that will be used. 
%
% In the working directory, make a subdirectory bak/matlab_owc/ and move
% the following files from the matlab_owc installation to there:
%
%	ow_config.txt
%	set_calseries.m
%
% These two files are now your own generic setup files for OWC, and you
% will edit them in (4-5). Then copies of them will be distributed into
% directories for each float (6) for more tailored settings as you
% analyse the float.

% 3) In your init_dmqc.m, edit your float list on top, as well as the
% information and paths below. Do not change the names of the objects,
% as these are tailored to be used by both DMQC-fun and OWC
% functions. The directory structure is built to keep the following
% types of files separated:
%
%	a) Matlab toolboxes
%	b) Data and figures produced by the toolboxes
%	c) Downloaded reference data
%	d) Downloaded float data
%	e) Your workspace with your own configurations and DMQC
%	   report material
%
% The philosphy is to be able to update toolboxes without losing
% settings or data (a), to store direct results but preferably be able
% to exclude them from any cloud sync or backup process which is both
% unneccessary and might also slow down your computer while processing
% (b), have the easily re-downloadable and reproduceable data in a
% temporary folder (c-d), but keeping your time consuming work and
% results where it is regularly backed up (e).  Ideas for this kind of
% placement can be found below on the lines marked 'EDIT THIS', where
% you wiill edit the paths according to your directory structure.

% 4) Edit the paths in your generic ow_config.txt. HISTORICAL_DIRECTORY,
% FLOAT_SOURCE_DIRECTORY, FLOAT_MAPPED_DIRECTORY, FLOAT_CALIB_DIRECTORY,
% FLOAT_PLOTS_DIRECTORY, CONFIG_DIRECTORY must be matched to
% owc_data_dir here in init_dmqc.m. Also, MAP_P_EXCLUDE should be the
% same in both files.

% 5) Edit the mapping and calibration parameters in your generic
% ow_config.txt and set_calseries.m. As far as possible and based on
% knowledge and advise from fellow operators in the region, set the
% scales etc. to something that should work in your region. Make
% comments and commented alternatives, so you later can keep track of
% what you are doing. The more alternatives and comments you put here,
% the easier it will be to fine-tune to the individual floats later.

% 6) Have INIT_DMQC build work environment for you by turning on the
% final part of this script and run the script. This will build a
% subfolder structure for all your floats and distribute copies of your
% generic config files from (5), as well as the LaTeX report files from
% the DMQC-fun toolbox, into each subfolder. You will now be able to
% configure mapping and calibration, make comments and write your
% summary/discussion, individually for each float. The other necessary
% pieces of information will be distrbuted by DMQC-functions as you use
% them. When adding new floats to your DMQC responsibility, you can run
% this again only selecting the new floats in the float selection part
% above. There is a safety mechanism to avoid overwriting existing float
% directories, but it is recommended that you turn off the building part
% once you have your float directories in place.

% 7) You are now ready to start doing core DMQC and salinity calibration
% on your floats. Head over to your copy of WORK_LOG for a description
% of the work flow.

% Finally some words about dependencies and files. The files you moved
% from the toolboxes could have been copied if you made absolutely sure
% which one is first in the path list. You may desire to keep the
% orignals in their place for backup/reference, but may run the risk of
% them being used instead of your copies. Careful use of commenting
% instead of deleting is recommended instead.  For the same reason
% plot_diagnostics_ow could be deleted from the matlab_owc toolbox, but
% it should be all right as long as the DMQC-fun toolbox has priority in
% your path list.  There are copies of ow_config.txt and set_calseries.m
% in the 'bak' folder of the DMQC-fun toolbox, but these are just
% examples of how to alter them. It is recommended to copy these files
% from the version of matlab_owc you have downloaded instead of here, to
% ensure compatibility.

% Now go carefully through the next sections and tailor to yourself.


% ------------- OPERATOR DATA ETC. (for report and D-files) ---------------					EDIT THIS!
dmqc_operator.name		= 'Jan Even Øie Nilsen'; 
dmqc_operator.parameter		= 'PRIMARY'; % I.e., CORE
dmqc_operator.orcid		= '0000-0003-2516-6106';
dmqc_operator.institution	= 'Institute of Marine Research (IMR)';
dmqc_operator.address		= 'Strandgaten 196, Bergen, Norway';
% Fill in more operators to the structure if needed, e.g.:
% dmqc_operator(2).name		= 'Joe Doe'; 
% dmqc_operator(2).parameter	= 'DOXY';
% dmqc_operator(2).orcid		= '2222-2222-2222-2222';
% dmqc_operator(2).institution	= 'University B';
% dmqc_operator(2).address	= 'Knowledge Street 123, Chippingford, Narnia';
user_manual_version = '3.41';
history_update = 'DMQC performed on CORE variables';


% -------------- ASSORTED PARAMETERS ---------------------------------------
MAP_P_EXCLUDE=400; % Keep as same as in ow_config.txt, but here for
                   % use in PREPARE_FLOATS. 
		   % [] In future versions this should be localised to
                   % float directory.


% -------------- PATHS (USE FULL PATHS AND DO NOT CHANGE THE NAME OF ANY OBJECTS): ---------------------------------------
% Many of these paths are used by both DMQC-fun and OWC-functions.
% MATLABOW directories:
owc_data_dir='/Users/a21627/Arkiv/data/matlab_owc/'; % Where you want OWC to put output.			EDIT THIS!
% In your copy of ow_config.txt,  HISTORICAL_DIRECTORY, FLOAT_SOURCE_DIRECTORY, FLOAT_MAPPED_DIRECTORY,
% FLOAT_CALIB_DIRECTORY, FLOAT_PLOTS_DIRECTORY, CONFIG_DIRECTORY must be edited according to your chosen owc_data_dir.
project_name='project_NorARGO/'; % Subdirectory to put in OWC's output directories.				EDIT THIS!
source_dir=[owc_data_dir,'float_source/'];			% Automatic. 
float_dirs=repmat({project_name},1,length(float_names));	% Automatic; necessary for ow_calibration.m.

% DOWNLOAD directories for float data:
download_parent_dir='/Users/a21627/Downloads/ARGO/'; % Where you want to download float data.			EDIT THIS!
download_dir=strcat(download_parent_dir,float_names,filesep);	% Automatic. 
rootdirin  = strcat(download_dir,'profiles/');			% Automatic. 
rootdirout = strcat(download_dir,'Dfiles/');			% Automatic; used by WRITE_D.

% REFERENCE DATA directories:
download_ref_data_dir='/Users/a21627/Downloads/DMQC'; % Where you want to download reference data.		EDIT THIS!
% Directories of reference data (versions):
struct2cell(dir(download_ref_data_dir)); ans(1,:); 
refdir=ans(contains(ans,'for_DMQC'));				% Automatic. 
% E.g., refdir={'CTD_for_DMQC_2021V01_1','CTD_for_DMQC_2021V01_7','ARGO_for_DMQC_2020V03'};
% Target directories (names of these subdirectories are given by OWC-toolbox):
tardir={[owc_data_dir,'climatology/historical_ctd'], ... 
	[owc_data_dir,'climatology/historical_bot'], ...
	[owc_data_dir,'climatology/argo_profiles']};		% Automatic. 

% REFERENCE DATA selection by WMO-squares:
% WMO-squares of chosen area to populate for DMQC-FUN and MATLAB_OWC.
% North Atlantic - Arctic sector:										EDIT THIS FOR YOUR REGION
my_WMOs = [1803 1804 1805	                         7802 7801 7800 1800 1801 1802 ...
	   1703 1704 1705	7707 7706 7705           7702 7701 7700 1700 1701 1702 ...
	   1603 1604 1605	     7606 7605 7604 7603 7602 7601 7600 1600 1601 1602 ...
			             7506 7505 7504 7503 7502 7501 7500 1500 1501 1502 ...
			             7406 7405 7404 7403 7402 7401 7400  ];
% After running LOAD_REFERENCEDATA you can check the map in matlab_owc's 'constants' directory.

% Your other relevant Argo directories:
my_argo_toolbox_dir = '/Users/a21627/matlab/toolbox/Argo/'; % Where you put all the Argo relevant toolboxes
							    % (i.e. DMQC-fun and MATLAB_OWC).			EDIT THIS!
my_working_dir      = '/Users/a21627/arbeid/obs/NorARGO/';  % Where you want to work with your Argo DMQC.	EDIT THIS!
my_backup_dir       = [my_working_dir,'bak/matlab_owc/'];   % Where you put your copies of ow_config_txt and 
							    % set_calseries.m.					EDIT THIS IF NECESSARY. 

% Make and (re)order paths in the necessary order (automatic, do not edit):
addpath([my_argo_toolbox_dir,'matlab_owc']);			% The OWC-toolbox' location
addpath([my_argo_toolbox_dir,'DMQC-fun']);			% This toolbox' location
addpath(my_working_dir);					% In order for your local version of INIT_DMQC to be used. 

% File names for PREPARE_FLOATS (automatic, do not edit):
infiles   = strcat(download_dir,filesep,float_names,{'_prof.nc'}); 
proffiles = strcat(download_dir,filesep,float_names,{'_prof.nc'}); 
techfiles = strcat(download_dir,filesep,float_names,{'_tech.nc'}); 
metafiles = strcat(download_dir,filesep,float_names,{'_meta.nc'}); 
trajfiles = strcat(download_dir,filesep,float_names,{'_Rtraj.nc'}); 
outfiles  = strcat(source_dir,float_dirs,float_names,'.mat');	% PREPARE_FLOATS uses this base filename also for plots.
								% (outfiles is also used by WRITE_D.)

% Directories for WRITE_D etc. (automatic, do not edit):
float_dir  = [owc_data_dir,'float_source',filesep,project_name]; % Also for DOWNLOAD_FLOATS and RUN_OW_CALIBRATION.
mapped_dir = [owc_data_dir,'float_mapped',filesep,project_name]; % Added for completeness.
calib_dir  = [owc_data_dir,'float_calib',filesep,project_name];
plot_dir   = [owc_data_dir,'float_plots',filesep,project_name]; % Also for DOWNLOAD_FLOATS (and the report)


% -------------- THE FIRST BUILD OF YOUR WORKING DIRECTORIES: -------------------- (automatic, do not edit)
if false      % (true = on, false = off)
	      % SWITCH ON _ONLY_ FOR THE _VERY_, _VERY_ FIRST TIME
	      % when you are building default workspaces for your set
              % of floats, OR when adding new floats to your set. In the
              % latter case, select only new floats in the float
              % selection section above. 
  mkdir([my_working_dir,'DMQC']);
  for I=1:length(float_names)
    [my_working_dir,'DMQC',filesep,float_names{I}];
    if ~exist(ans,'dir')
      disp(['Creating directory ',ans]);
      mkdir(ans);
      copyfile([my_backup_dir,'ow_config.txt'],  [my_working_dir,'DMQC',filesep,float_names{I},filesep]);
      copyfile([my_backup_dir,'set_calseries.m'],[my_working_dir,'DMQC',filesep,float_names{I},filesep]);
      copyfile([my_argo_toolbox_dir,'DMQC-fun/tex/DMQCreport.tex'],	[my_working_dir,'DMQC',filesep,float_names{I},filesep]);
      copyfile([my_argo_toolbox_dir,'DMQC-fun/tex/notes.tex'],		[my_working_dir,'DMQC',filesep,float_names{I},filesep]);
      copyfile([my_argo_toolbox_dir,'DMQC-fun/tex/discussion.tex'],	[my_working_dir,'DMQC',filesep,float_names{I},filesep]);
      copyfile([my_argo_toolbox_dir,'DMQC-fun/text/supplementary.tex'],	[my_working_dir,'DMQC',filesep,float_names{I},filesep]);
    end
  end
end
