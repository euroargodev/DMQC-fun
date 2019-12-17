# DMQC-fun
Extra scripts and functions to perform Argo Floats DMQC  

Why?
Between 2nd and 6th of December Even (Norway) and Gosia (Poland) visit Birgit and Ingrid in BSH (Germany) to learn the DMQC procedures. During the visit, we came up with a "wish list", where we listed the functionalities that we would like to add to the main OWC toolbox or that we need for the preparation of the D-files.

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

% Here is an overview of the files Even has made for a start. These
% files are concentrated around the pre-calibration stage. There are
% more documentation in the scripts themeselves.
%
% Note that some of them will need small functions from other toolboxes,
% in particular EVENMAT which can be downloaded at
% https://github.com/evenrev1/evenmat .
%
% init_dmqc	Name files to load/analyse, and set all paths for the 
%		DMQC. Download NetCDF file from coriolis server is
%		done here
% load referencedata	A script to ingest, quick-check and update
%		list of reference data for Argo DMQC.
% load_floats	Build the mat-file for the float, while checking for 
%		pressure and density, and existing QC-flags. 
% plot_profiles	Used by load_floats to make plots if any instabilities or
%		non-monotonic increasing pressure found.
% plot_diagnostics_ow	This is copied from the OWC toolbox and
%		modidified here, because we felt it needed some
%		improvements. It will be used instead of the
%		original/downloaded one if you run ow_calibration
%		from this local directory or make sure the path to
%		our version is first in Matlab's path list. 
%
% The following files are files in the OWC toolbox. I just mention these
% for completeness, as editing them by setting parameters in them are
% part of the recipe for calibration:
%
% ow_config.txt	Set the OWC paths and mapping parameters. Remember to
%		delete map- and calibration files for the float in
%		question. 
% set_calseries	Set the calibration parameters inside this
%		script. Remember to delete map- and calibration files
%		for the float in question. 
% ow_calibration Run the (mapping) and calibration.  
