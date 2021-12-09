% DMQC-fun 
% by J. Even Ø. Nilsen, Ingrid Angel, Birgit Klein, and Kjell Arne Mork.
% v0.9.3, jan.even.oeie.nilsen@hi.no.
%
% DMQC-fun is a comprehensive toolbox for performing DMQC on and
% salinity calibration of core data from Argo floats. The toolbox
% provides a system for semi-automated work-flow through the stages of
% downloading reference and float data, general DMQC and preparation of
% float data for MATLAB_OWC, running of MATLAB_WC, and production of
% updated D-files. Graphs, metadata and information from these stages
% are automaticly integrated into reports for each float.
%
% The toolboxes to install are:
% 
%	- DMQC-fun - https://github.com/imab4bsh/DMQC-fun.git
%	- matlab_owc -  https://github.com/ArgoDMQC/matlabow.git
%	- evenmat - https://github.com/evenrev1/evenmat.git 
%	- m_map - https://www.eoas.ubc.ca/~rich/map.html
%	- GSW - http://www.teos-10.org/pubs/gsw/html/
%
% DMQC-fun also includes a LaTeX report template and the functions and
% scripts produce snippets of content linked into that template. You
% will need a working version of LaTeX, and for instance dviPDFm to make
% PDF reports. However, this is not necessary for DMQC-fun to be
% useful. (You may even be able to link the produced figures and text
% parts into some other word processor of choice.)
%
% The files and directories in this toolbox are as follows:
%
% Contents	The Matlab help text for the toolbox.
%
% README.md	The Github version of Contents.m.
%
% init_dmqc	Both init-file and setup. Contains all installation
%		and setup instructions, as well as list of your
%		floats and parameters that all functions use. Read
%		that first! 
%
% work_log.txt	Explains the whole workflow of DMQC and the use of this
%		toolbox, and a good place to log your overall
%		progress. Understand more there! 
%
% load_referencedata	
%		A script to ingest, quick-check and update the list
%		of reference data for Argo DMQC and OWC. 
%
% download_floats	
%		This script downloads Argo float NetCDF-files from
%		the Coriolis server, as well as the altimetry
%		comparison and current greylist. 
%
% prepare_floats	
%		Does general DMQC, and builds the mat-files for OWC. 
%
% plot_profiles	Used by PREPARE_FLOATS to make plots if any instabilities or
%		non-monotonic increasing pressure is found.
%
% inpolygon_referencedata	
%		A function used by PREPARE_FLOATS to find reference
%		data inside a lon/lat polygon.
%
% operator_CPcor_new	
%		A tool for finding operator CPcor to correct for
%		pressure dependent conductivity bias, using near
%		float deployment ship CTD data.
%
% argo_sbe_CPcor	
%		A function used by OPERATOR_CPCOR_NEW to recaclulate
%		salinity with new CPcor.
%
% run_ow_calibration	
%		Runs OWC on all selected floats.
%
% write_D	Produces the D-files to deliver.
%
% plot_diagnostics_ow	
%		This function is copied from the OWC toolbox and
%		modified here because we felt it needed some
%		improvements. It will be used instead of the original
%		if you make sure the path to our version is first in 
%		Matlab's path list.
%
% ./tex/	The LaTeX part of the toolbox. You can ignore this
%		directory. It will be automaticly distributed when
%		you do your initial setup (see INIT_DMQC). And
%		there's a sample-DMQC-report.pdf here.
%
% ./lib/	The good old oceans and seawater toolboxes are
%		provided here since they are not easily found
%		anymore. Some other necessary tools are here as
%		well. Unzip and add paths as needed.
% 
% ./old/	Ignore this directory.
%
% The following files are only examples of how you can edit the two
% files from the OWC toolbox that you will copy to your working
% directory. But you can safely ignore this directory, and DO NOT add
% the path to it! See instructions in INIT_DMQC.
%
% ./bak/ow_config.txt	
%		Example of how to set the OWC paths and mapping
%		parameters.  
%
% ./bak/set_calseries.m	
%		Example of how to set the OWC calibration parameters.
%
% Feedback is most welcome to jan.even.oeie.nilsen@hi.no.
