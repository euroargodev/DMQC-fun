<<<<<<< HEAD
% DMQC-fun v0.9.1; 
% By J. Even Ø. Nilsen, Ingrid Angel, Birgit Klein, and Kjell Arne Mork.
%
% DMQC-fun is a comprehensive toolbox for performing DMQC on and
% salinity calibration of core data from Argo floats. The toolbox
% provides a system for semi-automated work-flow through the stages of
% downloading reference and float data, general DMQC and preparation of
% float data for MATLAB_OWC, and production of updated D-files. Graphs,
% metadata and information from these stages are automaticly integrated
% into reports for each float.
%
% The necessary toolbaoxes are:
% 
%	DMQC-fun - https://github.com/imab4bsh/DMQC-fun.git
%	matlab_owc -  https://github.com/ArgoDMQC/matlabow.git
%	evenmat - https://github.com/evenrev1/evenmat.git
%
% LaTeX: DMQC-fun includes a LaTeX report template and the Matlab
% scripts produce snippets of content linked into that template. You
% will need a working version of LaTeX, and for instance dvipdfm to make
% PDF. However, this is not necessary for DMQC-fun to be useful. (You
% may even be able to link the produced figures and text parts into some
% other word processor of choice.)
%
% The files in this toolbox are as follows:
%
% DMQC-fun	This Matlab help text for the toolbox.
%
% README.md	This file's twin on Github.
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
%		of reference data for Argo DMQC. 
%
% download_floats	
%		This script downloads float Argo NetCDF-files from
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
%		you do your initial setup (see INIT_DMQC).
%
% ./old/	Ignore this directory.
% 
% ./padconcatenation/
%		A copy of a useful auxillary toolbox (add to path list).
%
% The following files are only examples of how you can edit the two
% files from the OWC toolbox that you will copy to your working
% directory. But you can safely ignore this directory, and DO NOT add
% the path to it! See instructions in INIT_DMQC.
%
% ./bak/ow_config.txt	
%		Set the OWC paths and mapping parameters. 
%
% ./bak/set_calseries.m	
%		Set the OWC calibration parameters inside this
%		script.
%
% - eof -
=======
# DMQC-fun

NOTE: There is a possible major revision going on. Please hold ... -Even

Extra scripts and functions to perform Argo Floats DMQC  

Why?
Between 2nd and 6th of December Even (Norway) and Gosia (Poland) visit Birgit and Ingrid in BSH (Germany) to learn the DMQC procedures. During the visit, we came up with a "wish list", where we listed the functionalities that we would like to add to the main OWC toolbox or that we need for the preparation of the D-files.

See the matlab DMQC-fun.m documentation for matlab (help DMQC-fun) for a list of functions and more guidelines for installment and use. The scripts themselves (should) have the most updated instructions in for of good coomenting.

The List

1.	Script to write the D files (skeleton with annotation was send by Birgit to all of us)

2.	Plots
General: Plot only a certain section of the float (some cycles)
a.	Plot 1: Float trajectory (cycle color coded) and in blue the data used for the mapping
•	Distinguish between CTD and ARGO profiles
•	(On demand) Similar plot that shows which data was used for a given cycle
•	Plot bathymetry contours
•	Set legend outside of axes.
•	Control strictly the layout of eps-files (e.g., by setting Paper* properties of the matlab figure).
b.	Plot 4: Theta-S
•	Set TS axis ranges the same as Plot 3.
•	Control layout of eps-files (since there is a legend)
•	Plot reference data in one colour, mostly to limit the legend.
c.	Brian King plots
•	Define the ranges so it comprises the theta levels used for all cycles
•	Plot also negative values
•	Add a colorbar-colorcode next to the profile number to make the connection with the trajectory file
d.	Plot with the correction magnitude-piecewise regression
•	Add the 0.01 -0.01 green band
e.	Extra plot
On demand. For a given cycle, plot the the press vs. salinity, mapped salinity and color coded with year of the climatology (most recent)
f.	For the plot of calibrated data at different theta levels: select which levels to plot (now plots either 2 or 10)

3.	Ask Tierry Carval to add codes for operators-programs (Poland, Norway)

4.	Data massaging
a.	Apex floats: If pressure is corrected, then salinity has to be recalculated too (conductivity is calculated backwards from salinity and recalculated with new pressure)
b.	First check if the files are complying with the Real time (R mode rules), in particular: 
•	Pressure is not monotonically increasing;
•	Density inversions.
c.	Flag data (using scoop) to remove spikes-hooks

5.	Write a function that reads from R and D files instead of the *_prof.mat file, for making the mat file. (This is necessary if scoop works on the R and D files.)

6.	Creating the mat file for ow (aka further data massaging)

7.	Include functionality to avoid selecting the wrong data ( across basins )

8.	Write a function to assign grade to the entire profile according to the table in the manual (see skeleton for writing the d files)

9.	Write a function to put only the reference data you need for your region into the OWC climatology folders, including some checks and plots on it, as well as automatic update of your WMO box table in wmo_boxes-mat.

10.	Set up a github repository for our little group. INGRID
>>>>>>> 0a0b9bab027a992e9b016d6a67338278aaf432b4
