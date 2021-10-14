# DMQC-fun v0.9.1 
 by J. Even Ø. Nilsen, Ingrid Angel, Birgit Klein, and Kjell Arne Mork.

 DMQC-fun is a comprehensive toolbox for performing DMQC on and
 salinity calibration of core data from Argo floats. The toolbox
 provides a system for semi-automated work-flow through the stages of
 downloading reference and float data, general DMQC and preparation of
 float data for MATLAB_OWC, and production of updated D-files. Graphs,
 metadata and information from these stages are automaticly integrated
 into reports for each float.

 The necessary toolboxes are:
 
 - DMQC-fun - [https://github.com/imab4bsh/DMQC-fun.git](https://github.com/imab4bsh/DMQC-fun.git)
 - matlab_owc -  [https://github.com/ArgoDMQC/matlabow.git](https://github.com/ArgoDMQC/matlabow.git)
 - evenmat - [https://github.com/evenrev1/evenmat.git](https://github.com/evenrev1/evenmat.git) 

 DMQC-fun also includes a LaTeX report template and the functions and
 scripts produce snippets of content linked into that template. You
 will need a working version of LaTeX, and for instance dviPDFm to make
 PDF reports. However, this is not necessary for DMQC-fun to be
 useful. (You may even be able to link the produced figures and text
 parts into some other word processor of choice.)

 The files in this toolbox are as follows:

**Contents**	This Matlab help text for the toolbox.

**README.md**	This file's twin on Github.

**init_dmqc**	Both init-file and setup. Contains all installation
		and setup instructions, as well as list of your
		floats and parameters that all functions use. Read
		that first! 

**work_log.txt**	Explains the whole workflow of DMQC and the use of this
		toolbox, and a good place to log your overall
		progress. Understand more there! 

**load_referencedata**	
		A script to ingest, quick-check and update the list
		of reference data for Argo DMQC. 

**download_floats**	
		This script downloads float Argo NetCDF-files from
		the Coriolis server, as well as the altimetry
		comparison and current greylist. 

**prepare_floats**	
		Does general DMQC, and builds the mat-files for OWC. 

**plot_profiles**	Used by PREPARE_FLOATS to make plots if any instabilities or
		non-monotonic increasing pressure is found.

**inpolygon_referencedata**	
		A function used by PREPARE_FLOATS to find reference
		data inside a lon/lat polygon.

**run_ow_calibration**	
		Runs OWC on all selected floats.

**write_D**	Produces the D-files to deliver.

**plot_diagnostics_ow**	
		This function is copied from the OWC toolbox and
		modified here because we felt it needed some
		improvements. It will be used instead of the original
		if you make sure the path to our version is first in 
		Matlab's path list.

**./tex/**	The LaTeX part of the toolbox. You can ignore this
		directory. It will be automaticly distributed when
		you do your initial setup (see INIT_DMQC).

**./old/**	Ignore this directory.
 
**./padconcatenation/**	
		A copy of a useful auxillary toolbox (add to path list).

 The following files are only examples of how you can edit the two
 files from the OWC toolbox that you will copy to your working
 directory. But you can safely ignore this directory, and DO NOT add
 the path to it! See instructions in INIT_DMQC.

**./bak/ow_config.txt**	
		Set the OWC paths and mapping parameters. 

**./bak/set_calseries.m**	
		Set the OWC calibration parameters inside this
		script.

