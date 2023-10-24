% DOWNLOAD_FLOATS downloads float Argo NetCDF-files from the Coriolis
% server, as well as the altimetry comparison and current greylist.
%
% DMQC-fun v0.9.
% J. Even Ø. Nilsen, Ingrid M. Angel-Benavides, Birgit Klein, Malgorzata Merchel, and Kjell Arne Mork.
% Last updated: Wed Oct 18 11:45:42 2023 by jan.even.oeie.nilsen@hi.no

% You set which floats to operate on in INIT_DMQC! 
%
% No editing in this file (unless they change things on the ifremer-server).
%
% --------------------------------------------------------------------------

clear all; close all
init_dmqc; % float_names, paths, URLs, and filenames, etc.

for I=1:length(download_dir)
  if exist(download_dir{I},'dir'), rmdir(download_dir{I},'s'); end
  % Download the aggregated files:
  mkdir(download_dir{I}); 
  % FTP(host, username, password)
%fob=ftp(float_main_download_site);
  system(['lftp -e ''lcd ',download_dir{I},' ; get ',float_names{I},'_prof.nc ; bye'' ',float_main_download_site,float_names{I},'/']);
  system(['lftp -e ''lcd ',download_dir{I},' ; get ',float_names{I},'_tech.nc ; bye'' ',float_main_download_site,float_names{I},'/']);
  system(['lftp -e ''lcd ',download_dir{I},' ; get ',float_names{I},'_meta.nc ; bye'' ',float_main_download_site,float_names{I},'/']);
  system(['lftp -e ''lcd ',download_dir{I},' ; get ',float_names{I},'_Rtraj.nc ; bye'' ',float_main_download_site,float_names{I},'/']);
  system(['lftp -e ''lcd ',download_dir{I},' ; get ',float_names{I},'.png ; bye'' ftp.ifremer.fr/ifremer/argo/etc/argo-ast9-item13-AltimeterComparison/figures/']);
  % Download profiles from the R-files:
  mkdir(rootdirin{I}); 
  %system(['lftp -e ''lcd ',rootdirin{I},' ; mget D',float_names{I},'_*.nc ; bye'' ',float_profile_download_site,float_names{I},'/profiles/']);
  system(['lftp -e ''lcd ',rootdirin{I},' ; mget R',float_names{I},'_*.nc ; bye'' ',float_profile_download_site,float_names{I},'/profiles/']);
end

[download_parent_dir,filesep,'coriolis_greylist.csv']; if exist(ans,'file'), delete(ans); end
system(['lftp -e ''lcd ',download_parent_dir,' ; get coriolis_greylist.csv ; bye'' ftp.ifremer.fr/ifremer/argo/etc/greylist/']);


% WHAT THIS SCRIPT DOES, and what you might have to do yourself if
% your system does not allow to do this kind of batch download from
% command line. For each float do the following:
%
% 1) Make a local directory under download_parent_dir with your float's name as its name
% 2) Open the ftp site name in float_main_download_site in an FTP-client (anonymous login)
% 3) Go down the full path in float_main_download_site to find the directory with your float's name
% 3) Get and put the following files in the local folder (from 1) ending with:
%	_prof.nc
%	_tech.nc
%	_meta.nc
%	_Rtraj.nc
% 4) Find the png-figure in /ifremer/argo/etc/argo-ast9-item13-AltimeterComparison/figures/ and
%    put it in the same folder.
% 5) Then make a local subfolder called 'profiles' inside the same folder.
% 6) Go to the path in float_profile_download_site and get all the files starting with R (i.e. R*.nc)
% 7) Also get coriolis_greylist.csv from the directory /ifremer/argo/etc/greylist/



% R-files must be downloaded at the same time, as done here, or else there might be mismatch between number of cal parameters and R-files.


