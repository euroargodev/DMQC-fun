% DOWNLOAD_FLOATS downloads float Argo NetCDF-files from the Coriolis
% server, as well as the altimetry comparison and current greylist.
% by J. Even Ø. Nilsen, Ingrid Angel, Birgit Klein, and Kjell Arne Mork.
% DMQC-fun v0.9.3, jan.even.oeie.nilsen@hi.no.
% 
% You set which floats to operate on in INIT_DMQC!

clear all; close all
init_dmqc; % Paths and filenames, etc.

for I=1:length(download_dir)
  if exist(download_dir{I},'dir'), rmdir(download_dir{I},'s'); end
  mkdir(download_dir{I}); 
  system(['lftp -e ''lcd ',download_dir{I},' ; get ',float_names{I},'_prof.nc ; bye'' ftp.ifremer.fr/ifremer/argo/dac/coriolis/',float_names{I},'/']);
  system(['lftp -e ''lcd ',download_dir{I},' ; get ',float_names{I},'_tech.nc ; bye'' ftp.ifremer.fr/ifremer/argo/dac/coriolis/',float_names{I},'/']);
  system(['lftp -e ''lcd ',download_dir{I},' ; get ',float_names{I},'_meta.nc ; bye'' ftp.ifremer.fr/ifremer/argo/dac/coriolis/',float_names{I},'/']);
  system(['lftp -e ''lcd ',download_dir{I},' ; get ',float_names{I},'_Rtraj.nc ; bye'' ftp.ifremer.fr/ifremer/argo/dac/coriolis/',float_names{I},'/']);
  system(['lftp -e ''lcd ',download_dir{I},' ; get ',float_names{I},'.png ; bye'' ftp.ifremer.fr/ifremer/argo/etc/argo-ast9-item13-AltimeterComparison/figures/']);
  %crop([download_dir{I},float_names{I},'.png']);
  % pres_adj
  %system(['lftp -e ''lcd ',download_dir{I},' ; get ',float_names{I},'_BRtraj.nc ; bye'' ftp.ifremer.fr/ifremer/argo/dac/coriolis/',float_names{I},'/']);
  %system(['lftp -e ''lcd ',download_dir{I},' ; get ',float_names{I},'_Sprof.nc ; bye'' ftp.ifremer.fr/ifremer/argo/dac/coriolis/',float_names{I},'/']);
  %
  % From May 2021, download profiles from the R-files:
  mkdir(rootdirin{I}); 
  system(['lftp -e ''lcd ',rootdirin{I},' ; mget D',float_names{I},'_*.nc ; bye'' ftp.ifremer.fr/ifremer/argo/dac/coriolis/',float_names{I},'/profiles/']);
  system(['lftp -e ''lcd ',rootdirin{I},' ; mget R',float_names{I},'_*.nc ; bye'' ftp.ifremer.fr/ifremer/argo/dac/coriolis/',float_names{I},'/profiles/']);
end

[download_parent_dir,filesep,'coriolis_greylist.csv']; if exist(ans,'file'), delete(ans); end
system(['lftp -e ''lcd ',download_parent_dir,' ; get coriolis_greylist.csv ; bye'' ftp.ifremer.fr/ifremer/argo/etc/greylist/']);


% R-files must be downloaded at the same time, as done here, or else there might be mismatch between number of cal parameters and R-files.


