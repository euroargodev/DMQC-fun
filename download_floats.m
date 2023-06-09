% DOWNLOAD_FLOATS downloads float Argo NetCDF-files from the Coriolis
% server, as well as the altimetry comparison and current greylist.
%
% DMQC-fun v0.9.
% J. Even Ø. Nilsen, Ingrid M. Angel-Benavides, Birgit Klein, Malgorzata Merchel, and Kjell Arne Mork.
% Last updated: Wed May 24 13:40:44 2023 by jan.even.oeie.nilsen@hi.no

% You set which floats to operate on in INIT_DMQC! 
%
% No editing in this file (unless they change things on the ifremer-server).
%
% --------------------------------------------------------------------------

clear all; close all
init_dmqc; % Paths and filenames, etc.

for I=1:length(download_dir)
  if exist(download_dir{I},'dir'), rmdir(download_dir{I},'s'); end
  % Download the aggregated files:
  mkdir(download_dir{I}); 
  system(['lftp -e ''lcd ',download_dir{I},' ; get ',float_names{I},'_prof.nc ; bye'' ',float_main_download_site,float_names{I},'/']);
  system(['lftp -e ''lcd ',download_dir{I},' ; get ',float_names{I},'_tech.nc ; bye'' ',float_main_download_site,float_names{I},'/']);
  system(['lftp -e ''lcd ',download_dir{I},' ; get ',float_names{I},'_meta.nc ; bye'' ',float_main_download_site,float_names{I},'/']);
  system(['lftp -e ''lcd ',download_dir{I},' ; get ',float_names{I},'_Rtraj.nc ; bye'' ',float_main_download_site,float_names{I},'/']);
  system(['lftp -e ''lcd ',download_dir{I},' ; get ',float_names{I},'.png ; bye'' ftp.ifremer.fr/ifremer/argo/etc/argo-ast9-item13-AltimeterComparison/figures/']);
  % Download profiles from the R-files:
  mkdir(rootdirin{I}); 
  system(['lftp -e ''lcd ',rootdirin{I},' ; mget D',float_names{I},'_*.nc ; bye'' ',float_profile_download_site,float_names{I},'/profiles/']);
  system(['lftp -e ''lcd ',rootdirin{I},' ; mget R',float_names{I},'_*.nc ; bye'' ',float_profile_download_site,float_names{I},'/profiles/']);
end

[download_parent_dir,filesep,'coriolis_greylist.csv']; if exist(ans,'file'), delete(ans); end
system(['lftp -e ''lcd ',download_parent_dir,' ; get coriolis_greylist.csv ; bye'' ftp.ifremer.fr/ifremer/argo/etc/greylist/']);


% R-files must be downloaded at the same time, as done here, or else there might be mismatch between number of cal parameters and R-files.


