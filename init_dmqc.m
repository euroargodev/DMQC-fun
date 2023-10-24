% INIT_DMQC.M is init-file and setup instructions.  
% See Contents.m for overview; work_log.txt for workflow description.
%
%%%%%%%%%% Read further down below for instructions on how to set up your system! %%%%%%%%%%%%%%%%  
%%%%%%%%%% The first part is where you control everything! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DMQC-fun v0.9.
% J. Even Ø. Nilsen, Ingrid M. Angel-Benavides, Birgit Klein, Malgorzata Merchel, and Kjell Arne Mork.
% Last updated: Tue Oct 24 12:16:27 2023 by jan.even.oeie.nilsen@hi.no

% ------------ FLOAT SETTINGS: ----------------------------------------
TCL=complex(0,1); % Do not touch!
% List all your floats here:
	     
%             1         2         3         4         5         6         7         8         9         10      
float_names={'6903549','6903550','6903551','6903552','6903553','6903554','6903555','6903556','6903557','6903558',... % 0
	     '6903559','6903560','6903561','6903562','6903563','6903564','6903565','6903566','6903567','6903568',... % 10
	     '6903569','6903570','6903571','6903572','6903573','6903574','6903575','6903577','6903578','6903579',... % 20
	     '6903580','6903581','6903582','6903583','6903584','6903585','6903586','6903587','6903588','6903589',... % 30
	     '6903576','6903590','6903591','6903592','3902462','3902463','3902464','7901006','7901007','3902465',... % 40
	     '1902579','2903771','6904242'									};   % 50
%             1         2         3         4         5         6         7         8         9         10      
cal_action ={ [0]     , [0]     , [1]     , [0]     , [0]     , [0]     , [0]     , [0 1 4] , [0 1 4] , [0]     ,... % 0 
	      [0]     , [0]     , [0 1]   , [0 1 4] , [0]     , [0]     , [0]     , [0]     , [0]     , [1]     ,... % 10
	      [0]     , [0]     , [0]     , [0]     , [0]     , [1]     , [0]     , [0]     , [0]     , [0 0]   ,... % 20 
	      [0]     , [0]     , [1]     , [0]     , [1]     , [1]     , [1]     , [0]     , [0]     , [0]     ,... % 30
	      [0]     , [1]     , [1]     , [0]     , [0]     , [0]     , [0]     , [0]     , [0]     , [0]     ,... % 40
	      [0]     , [0]     , [1]										};   % 50
%             1         2         3         4         5         6         7         8         9         10      
checked    ={ 230       147       124       52        113       152       70        100       124       127     ,... % 0 
	      0         0         136       132       0         412       369       0         176       89      ,... % 10
	      164       185       104       90        90        186       139       141       143       154     ,... % 20 
	      67        201       5         0         0         0         0         0         162       0       ,... % 30
	      47        43        40        47        22        23        23        45        54        23      ,... % 40
              23        21        64										};   % 50
%             1         2         3         4         5         6         7         8         9         10      
Dchecked   ={ 0         0         0         0         0         0         0         1         1         1       ,... % 0 
	      0         0         1         1         0         1         1         0         0         0       ,... % 10
              0         0         1         1         1         87        139       141       143       154     ,... % 20 
	      0         1         0         0         0         0         0         0         1         0       ,... % 30
	      47        43        40        47        1         1         1         45        21        1       ,... % 40
              0         1         1										};   % 50
%             1         2         3         4         5         6         7         8         9         10      
Nclu       ={ 2         2         2         1         1         3         1         2*TCL     2         3*TCL   ,... % 0 
	      2         2         3*TCL     3*TCL     2         2         2         1         2         1       ,... % 10
	      2         2         2*TCL     2         1         1         1         2         1         {75 135},... % 20 
	      1         {70 135}  2         2         1         1         2         1         {40 120}  1	,... % 30     
	      1         2*TCL     1         1         1         1         1         1         1         1       ,... % 40
              1         1         1										};   % 50
%             1         2         3         4         5         6         7         8         9         10      

% ------------ FLOAT SELECTION: ----------------------------------------
%
% You can select which floats to analyse by just giving specific index
% numbers and all the onjects will be reduced to your
% selection. Helping-numbers for finding indices are given in the
% commented lines and end of lines above.  
% Examples:
1:length(float_names);	% All floats
[11:18 32:40 51:53]  	% Core (Arvor)
% [1:3 26:27 41:44]	% BGC (PROVOR CTS4)
[28:30 48:49]		% Bio (PROVOR CTS4)
% [4:7 19:22]		% Bio (APEX AFP11) 
% [8:10 23:25 31 45:47]	% Dyp (Arvor Deep) 
% [50]			% Core + DO (Arvor)
% %[34]			% √ Lost in the Arctic?
% setdiff(ans,[34])	% Exclude lost in the Arctic?
% [8 9 13 14];		% PSAL greylisted

ans(1) % To just look at the first of the selected group

% DO NOT CHANGE:
float_names=float_names(ans); Nclu=Nclu(ans); cal_action=cal_action(ans); checked=checked(ans);  Dchecked=Dchecked(ans); 
float_names % Display which floats are about to be processed

% ------------ DIRECTION SELECTION: ------------------------------------
%
% Indicate by commenting lines here whether you are doing DMQC on
direction = 'A'; % ascending profiles (Do these first!) or
%direction = 'D'; % descending profiles (do these separately after OWC).
% All the direction related changes are done in each script separately. 

% ------------ REFERENCE DATA SELECTION: -------------------------------
%
% Indicate by commenting lines here whether you are going to use:
itardir=[1 2 3]; % All reference data (DEFAULT)
%itardir=[1   3]; % only CTD and Argo
%itardir=[1    ]; % only CTD
%itardir=[    3]; % only Argo
% This will affect plots made by DMQC-fun (sets tardir), as well as
% OWC (wmo_boxes.mat will be replaced).


% ------------ EXPLANATION OF VARIABLES ABOVE: -------------------------
%
% float_names is necessary for ow_calibration.m. DO NOT change the
%	name of this variable!
% cal_action is your calibration decision. Set this initially to [0]
%	for every float.
%	case 0	Data are good; no adjustment to be applied.
%	case 1	Data show sensor drift or offset; adjustment to be applied.
%	case 4	% Data are bad and unadjustable.
% checked is the number of ascending profiles already manually checked
%	(you will be prompted to update this number in the
%	process). Set this initially to 0 for every float. 
% Dchecked is the number of DESCENDING profiles already manually
%	checked (you will be prompted to update this number in the
%	process). Set this initially to 0 for every float. 
% Nclu is specification of the clustering of profiles to use in
%	the visual comparison with reference data in PREPARE_FLOATS:
%	scalar - the number of clusters. Set to 2 or more if float
%		has traversed several basins or regions. 
%       scalar multiplied with TCL - means that time (i.e., cycle
%		number) will affect the clustering.
%	cell - a series of specific profile (!) numbers that will
%		separate the clusters, in order to control the
%		clusters manually.
%	char - any character means that the clustering will follow
%		the grouping defined by calseries as set in
%		set_calseries.m for the float. This advanced option
%		is only to be used in combination with deliberate 
%		grouping of the calseries. For descending profiles,
%		for which the calseries do not necessarily match,
%		there is possibility to indicate number of clusters
%		by choice of character (A=1, B=2, etc.).
% direction is descrived above.
% itardir is descrived above.

% ------------- OVERVIEW OF THIS INIT FILE: --------------------------
%
% This file will be used by all the main DMQC-fun scripts, as well as
% OWC functions, when using the DMQC-fun setup. It has the following
% sections:
% - A list of all your floats and corresponding parameters you will
%   change during DMQC (above for convenience). 
% - Float selection part for which floats to DMQC (above).
% - Direction selection for working on ascending or descending
%   profiles (above).  
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
% Argo related toolboxes (e.g. /Users/a21627/matlab/toolbox/Argo/):
% 
%	DMQC-fun - https://github.com/euroargodev/DMQC-fun.git
%	matlab_owc -  https://github.com/ArgoDMQC/matlabow.git
% 
% Also install in your general toolbox directory
% (e.g. /Users/a21627/matlab/toolbox/):
%
%	evenmat - https://github.com/evenrev1/evenmat.git
%	m_map - https://www.eoas.ubc.ca/~rich/map.html
%	GSW - http://www.teos-10.org/pubs/gsw/html/
%
% The oceans, seawater, and some other necessary toolboxes are provided
% herein. Unzip and add their paths if you don't have them installed
% already.
%
% As indicated, best practice of installing matlab toolboxes is to make
% a directory in your home folder called ~/matlab/toolbox/ and put all
% the toolbox folders there. Also, using DMQC-fun you should gather all
% your Argo toolboxes under a common directory ~/matlab/toolbox/Argo/.
%
% Then add paths to all the new toolboxes. Do not include subfolders at
% first, only do that later when error messages dictates which specific
% subfolders are needed.  Also, make sure paths to DMQC-fun comes above
% paths to the other toolboxes in the path list (i.e. addpath it last),
% so that our altered versions of some files (e.g. plot_diagnostics_ow)
% are used instead of the matlab_owc version. However, there is a
% fail-safe mechanism buildt into INIT_DMQC that adds and orders these
% paths upon first time use.
% 
% LaTeX: DMQC-fun includes a LaTeX report template and the Matlab
% scripts produce snippets of content linked into that template. You
% will need a working version of LaTeX, and for instance dvipdfm to make
% PDFs. However, this is not necessary for DMQC-fun to be working or
% even useful. (You may even be able to link the produced figures and
% text parts into some other word processor of choice.)

% 2) You need to copy some files from the toolboxes manually, as you do
% not want updating of toolboxes to erase your own settings. From
% DMQC-fun toolbox copy the following files to your working directory
% for Argo DMQC (e.g. ~/mywork/Argo/) (i.e. my_working_dir):
%
%	init_dmqc.m (this file) 
%	work_log.txt
%
% A path to your working directory will be added by INIT_DMQC itself to
% make sure this local INIT_DMQC is the one that will be used.
%
% In the working directory, make a subdirectory bak/matlab_owc/ and copy
% the following files from the matlab_owc installation to there:
%
%	ow_config.txt
%	set_calseries.m
%
% These two files are now your own generic setup files for OWC, and you
% will edit them (4-5). Then copies of them will be distributed into
% directories for each float (6) for more tailored settings as you
% analyse the float.

% 3) In your init_dmqc.m (your copy of this file), edit your float list
% on top, as well as the information and paths below. Do not change the
% names of the variables, as these are tailored to be used by both
% DMQC-fun and OWC functions. The directory structure is built to keep
% the following types of files separated:
%
%	a) Matlab toolboxes
%	b) Data and figures produced by the toolboxes
%	c) Downloaded reference data
%	d) Downloaded float data
%	e) Your workspace with your own configurations and DMQC
%	   report material
%
% Example structures, with your home folder as base:
% a) 
% ~/matlab/toolbox/			(all other toolboxes)
% ~/matlab/toolbox/Argo/		(my_argo_toolbox_dir; all Argo toolboxes)
% ~/matlab/toolbox/Argo/DMQC-fun/
% ~/matlab/toolbox/Argo/matlab_owc/
% b)
% ~/Arkiv/data/matlab_owc/		(owc_data_dir)
% c)
% ~/Downloads/DMQC/			(download_ref_data_dir)
% d) 
% ~/Downloads/ARGO/			(download_parent_dir)
% e)
% ~/mywork/Argo/			(my_working_dir)
% ~/mywork/Argo/DMQC/			(the individual float folders)
% ~/mywork/Argo/bak/			(any files you want to safekeep for reference)
% ~/mywork/Argo/bak/matlab_owc/		(your generic OWC-setup files; see (2))
% 
% The philosophy behind this separation is to be able to update
% toolboxes without losing settings or data (a), to store direct results
% but preferably be able to exclude them from any cloud sync or backup
% process or virus scanning which is unneccessary and might also
% slow down your computer while processing (b), have the easily
% re-downloadable and reproduceable data in a temporary folder (c-d),
% but keeping your time consuming work and results where it is regularly
% backed up (e).  Ideas for this kind of placement can be found below on
% the lines marked 'EDIT THIS', where you will edit the paths according
% to your own directory structure.

% 4) Edit the paths in your own copy of ow_config.txt (see point 2
% above). HISTORICAL_DIRECTORY, FLOAT_SOURCE_DIRECTORY,
% FLOAT_MAPPED_DIRECTORY, FLOAT_CALIB_DIRECTORY, FLOAT_PLOTS_DIRECTORY,
% CONFIG_DIRECTORY must be matched to owc_data_dir as set below here in
% your init_dmqc.m.

% 5) Edit the mapping and calibration parameters in your generic
% ow_config.txt and set_calseries.m. As far as possible and based on
% knowledge and advise from fellow operators in the region, set the
% scales etc. to something that should in general work in your
% region. Make comments and commented alternatives, so you later can
% keep track of what you are doing. The more alternatives and comments
% you put here, the easier it will be to fine-tune to the individual
% floats later.

% 6) Have INIT_DMQC build your work environment for you by setting float
% selection to include all floats and turning on the THE FIRST BUILD OF
% YOUR WORKING DIRECTORIES part at the end of this script, and run the
% script. This will build a subfolder structure for all your floats and
% distribute copies of your generic config files from (5), as well as
% the LaTeX report files from the DMQC-fun toolbox, into each subfolder,
% one for each float. You will now be able to configure mapping and
% calibration, make comments and write your summary/discussion,
% individually for each float. The other necessary pieces of information
% will be distrbuted by DMQC-functions as you use them. When adding new
% floats to your DMQC responsibility, you can create new folders for
% them by _carefully_ selecting _only_ these new floats in the float
% selection part above, turn on the building part, and run this
% script. There is a safety mechanism to avoid overwriting existing
% float directories, but it is recommended that you turn off the
% building part once you have your float directories in place.

% 7) You are now ready to start doing core DMQC and salinity calibration
% on your floats. Head over to your copy of WORK_LOG for a description
% of the work flow.

% Finally some words about dependencies and files. We copy files out of
% the toolboxes and into our own system, instead of moving, in order to
% keep the orignals in their place for backup/reference. However, then
% you must make absolutely sure which one is first in the path list. If
% the original toolboxes are listed above the DMQC-toolbox and your
% local system in the path list, the originals will be used
% instead. This must be avoided by careful attention to the path list.
% There are copies of ow_config.txt and set_calseries.m in the 'bak'
% folder of the DMQC-fun toolbox, but these are just examples of how to
% alter them. It is recommended to copy these files from the version of
% matlab_owc you have downloaded instead of from here, to ensure
% compatibility.

% Now go carefully through the next sections and tailor it to yourself.


% ------------- OPERATOR DATA ETC. (for report and D-files) ------------					EDIT THIS!
dmqc_operator.name		= 'Jan Even Øie Nilsen'; 
dmqc_operator.parameter		= 'PRIMARY'; % I.e., CORE
dmqc_operator.orcid		= '0000-0003-2516-6106';
dmqc_operator.institution	= 'Institute of Marine Research (IMR)';
dmqc_operator.address		= 'Strandgaten 196, Bergen, Norway';
% Fill in more operators to the structure if needed, e.g.:
% dmqc_operator(2).name		= 'Joe Doe'; 
% dmqc_operator(2).parameter	= 'DOXY';
% dmqc_operator(2).orcid	= '2222-2222-2222-2222';
% dmqc_operator(2).institution	= 'University B';
% dmqc_operator(2).address	= 'Knowledge Street 123, Chippingford, Narnia';
user_manual_version		= '3.41';
history_update			= 'DMQC performed on CORE variables';


% -------------- ASSORTED PARAMETERS -----------------------------------
% None at the moment.

% -------------- FTP SITES ----------------------------------------------					EDIT THIS FOR YOUR DAC!
float_main_download_site='ftp.ifremer.fr/ifremer/argo/dac/coriolis/'; % Where we normally download from, and where our D-files go when submitted.
%float_profile_download_site='ftp.ifremer.fr/ifremer/argo/dac/coriolis/'; % Location of 'profiles' directory (normally same place).
float_profile_download_site='ftp.ifremer.fr/ifremer/coriolis/argo/dac/coriolis/'; % EDAC
% EDAC is a special place if you need to access R files. These files are
% R-files that have format updated, but does not seem to have RTQC
% applied to them. Quoting Christine C.: "directly from our EDAC which
% provides the last version of the Argo float files" ; "the best way is
% to use directly the EDAC because you are sure to get the last version
% of the files.".

% -------------- PATHS (USE FULL PATHS AND DO NOT CHANGE THE NAME OF ANY OBJECTS): ---------------------------------------
% Many of these paths are used by both DMQC-fun and OWC-functions.
% MATLABOW directories:
owc_data_dir='/Users/a21627/Arkiv/data/matlab_owc/';		% Where you have set OWC to put output.		EDIT THIS!
% In your copy of ow_config.txt,  HISTORICAL_DIRECTORY, FLOAT_SOURCE_DIRECTORY, FLOAT_MAPPED_DIRECTORY,
% FLOAT_CALIB_DIRECTORY, FLOAT_PLOTS_DIRECTORY, CONFIG_DIRECTORY must be edited according to your chosen owc_data_dir.
project_name='project_NorARGO/'; % Subdirectory you have put in OWC's output directories (empty if not used).	EDIT THIS!
source_dir=[owc_data_dir,'float_source/'];			% Automatic. 
float_dirs=repmat({project_name},1,length(float_names));	% Automatic; necessary for ow_calibration.m.

% DOWNLOAD directories for float data:
download_parent_dir='/Users/a21627/Downloads/ARGO/'; % Where you want to download float data.			EDIT THIS!
download_dir=strcat(download_parent_dir,float_names,filesep);	% Automatic. 
rootdirin  = strcat(download_dir,'profiles/');			% Automatic. 
rootdirout = strcat(download_dir,'Dfiles/');			% Automatic; used by WRITE_D.

% REFERENCE DATA directories:
download_ref_data_dir='/Users/a21627/Downloads/DMQC'; % Where you will download reference data.			EDIT THIS!
% Directories of reference data (versions):
struct2cell(dir(download_ref_data_dir)); ans(1,:); 
refdir=ans(contains(ans,'_for_DMQC_'));				% Automatic. 
% E.g., refdir={'CTD_for_DMQC_2021V01_1','CTD_for_DMQC_2021V01_7','ARGO_for_DMQC_2020V03'};
% Target directories (names of these subdirectories are given by OWC-toolbox):
tardir={[owc_data_dir,'climatology/historical_ctd'], ... 
	[owc_data_dir,'climatology/historical_bot'], ...
	[owc_data_dir,'climatology/argo_profiles']};		% Automatic. 
tartyp={'ctd','bottle','argo'};					% Automatic. 

% REFERENCE DATA selection by WMO-squares:
% WMO-squares of chosen area to populate for DMQC-FUN and MATLAB_OWC.
% North Atlantic - Arctic sector:										EDIT THIS FOR YOUR REGION
my_WMOs = [                         7802 7801 7800 1800 1801 1802 1803 1804 1805 1806 1807 1808	...
	   7707 7706 7705           7702 7701 7700 1700 1701 1702 1703 1704 1705 1706 1707 1708	...
	        7606 7605 7604 7603 7602 7601 7600 1600 1601 1602 1603 1604			...
	        7506 7505 7504 7503 7502 7501 7500 1500 1501 1502				...
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
  [my_working_dir,'DMQC']; if ~exist(ans,'dir') mkdir(ans); end
  for I=1:length(float_names)
    [my_working_dir,'DMQC',filesep,float_names{I}];
    if ~exist(ans,'dir')
      disp(['Creating directory ',ans]);
      mkdir(ans);
      copyfile([my_backup_dir,'ow_config.txt'],[my_working_dir,'DMQC',filesep,float_names{I},filesep]);
      copyfile([my_backup_dir,'set_calseries.m'],[my_working_dir,'DMQC',filesep,float_names{I},filesep]);
      copyfile([my_argo_toolbox_dir,'DMQC-fun/tex/DMQCreport_float.tex'],[my_working_dir,'DMQC',filesep,float_names{I},filesep]);
      copyfile([my_argo_toolbox_dir,'DMQC-fun/tex/notes.tex'],[my_working_dir,'DMQC',filesep,float_names{I},filesep]);
      copyfile([my_argo_toolbox_dir,'DMQC-fun/tex/discussion.tex'],[my_working_dir,'DMQC',filesep,float_names{I},filesep]);
    end
  end
end

% --------------- FOR JUST DISTRIBUTING NEW VERSIONS OF LaTeX MASTER FILE: -------------------------------
if false      % (true = on, false = off)
  for I=1:length(float_names)
      copyfile([my_argo_toolbox_dir,'DMQC-fun/tex/DMQCreport_float.tex'],[my_working_dir,'DMQC',filesep,float_names{I},filesep]);
      % It is not recommended to distribute new notes.tex or
      % discussion.tex, as these will likely contain important information.
  end
end
