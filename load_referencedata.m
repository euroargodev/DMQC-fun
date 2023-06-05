% LOAD_REFERENCEDATA A script to ingest, quick-check and update list of
% reference data for Argo DMQC.
%
% DMQC-fun v0.9.
% J. Even Ø. Nilsen, Ingrid M. Angel-Benavides, Birgit Klein, Malgorzata Merchel, and Kjell Arne Mork.
% Last updated: Wed May 24 13:40:44 2023 by jan.even.oeie.nilsen@hi.no

% - Ingests downloaded reference data into MATLAB_OWC's climatology/ directory.
% - Creates three versions of the matrix of available metadata
%   (wmo_boxes.mat), one for each type in addition to the main file.
% - Produces map overview of the chosen WMO squares and contents in constants/)
% - Produces maps and TS, T, and S graphs of all profiles for each
%   ingested WMO-square alongside their mat-files. 
%
%%%%%%%%%%%% INITIAL STEPS: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First, in a shell, find, download and unpack DMQC reference data:
%	cd ~/Downloads/DMQC   (i.e., your download_ref_data_dir)
%	lftp -u <user,password> ftp.ifremer.fr/coriolis/
% Check for new datasets, get the tarbals if any.
% Close lftp session and unpack the tarballs inside your download_ref_data_dir.
%
% IMPORTANT: INIT_DMQC will build the object refdir according to the
% subdirectories present in the download_ref_data_dir, so make sure to
% delete old datasets and keep only the ones to be used!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% After running LOAD_REFERENCEDATA you can check the map in
% matlab_owc's 'constants' directory.
%
% No editing in this file!
%
% ------------------------------------------------------------------------------

clear all; close all;
init_dmqc;	% Get the necessary paths etc. as well as your
                % choice of WMO-squares.

% ----- Make an overview map of selected and ingested WMO-squares -----
figure(1); clf; set(gcf,'OuterPosition',get(0,'ScreenSize'));
wmolim=[]; for i=1:length(my_WMOs), wmolim(i,:)=wmosquare(my_WMOs(i)); end
mima(wmolim(:,3:4))+[-2 2];ans(2)=min(ans(2),89);
m_proj('Albers','lon',mima(wmolim(:,1:2))+[-2 2],'lat',ans);
m_grid; m_coast('color','k'); m_elev('contour','color',[.7 .7 .7]);
for i=1:length(my_WMOs)
  [lo,la]=erect(wmolim(i,:),'t'); m_line(lo,la);
  sort([wmolim(i,1:2)+[-4 4],wmolim(i,1:2)]);lom=ans(2:3);
  sort([wmolim(i,3:4)+[-1 1],wmolim(i,3:4)]);lam=ans(2:3);
  htx=m_text(mean(lom),mean(lam),int2str(my_WMOs(i))); 
  set(htx,'verticalalignment','middle','horizontalalignment','center','color','b');
end
% map will be printed after the content check below.


% ----- Transfer the reference data to the apropriate directories: -----
for j=1:length(refdir) % Loop reference-data directories
  if any(findstr('CTD',refdir{j}))
    for i=1:length(my_WMOs)
      %['cp ',download_ref_data_dir,filesep,refdir{j},filesep,'ctd_',int2str(my_WMOs(i)),'.mat ',tardir{1},filesep]
      msg=system(['cp ',download_ref_data_dir,filesep,refdir{j},filesep,'ctd_',int2str(my_WMOs(i)),'.mat ',tardir{1},filesep]);
    end
  elseif any(findstr('ARGO',refdir{j}))
    for i=1:length(my_WMOs)
      %['cp ',download_ref_data_dir,filesep,refdir{j},filesep,'argo_',int2str(my_WMOs(i)),'.mat ',tardir{3},filesep]
      msg=system(['cp ',download_ref_data_dir,filesep,refdir{j},filesep,'argo_',int2str(my_WMOs(i)),'.mat ',tardir{3},filesep]);
    end
  end
end


% ----- Update the wmo_boxes.mat based on actual contents of target directories: -----
%
% The 1st column in wmo_boxes.mat is the list of all WMO box
% numbers. The 2nd column denotes existence of CTD data. The 3rd column
% denotes existence of BOTTLE data. The 4th column denotes existence of
% Argo data. 0 = no data, or do not use. 1 = data exist, and use
% them. Edit the 2nd, 3rd and 4th columns after you have added the
% reference data.
load([owc_data_dir,filesep,'constants',filesep,'wmo_boxes']);
N=size(la_wmo_boxes,1);
% Search all target directories and update the 'flags':
la_wmo_boxes(:,2:4)=0;					% Reset
for j=1:length(tartyp)	% Loop the three data types
  d=dir(tardir{j}); filer=cellstr(char(d.name));	% List of files
  for i=1:N		% Loop all lines in the WMO-boxes list
    if any(contains(filer,[tartyp{j},'_',int2str(la_wmo_boxes(i,1)),'.mat']))
      la_wmo_boxes(i,1+j)=1;
      % Add text to the map for which data types exist in the WMOsq.
      wmosquare(la_wmo_boxes(i,1)); 
      la=linspace(ans(3),ans(4),10); la=la(1+j);
      htxex=m_text((ans(1)+ans(2))/2,la,tartyp{j});
      set(htxex,'HorizontalAlignment','center');
    end
  end
end
print(gcf,'-depsc',[owc_data_dir,filesep,'constants',filesep,'selected_wmos_map.eps']);
save([owc_data_dir,filesep,'constants',filesep,'wmo_boxes'],'la_wmo_boxes');

% Save alternative versions in order to selectively OWC with CTD or Argo only:
la_wmo_boxes_ctd=la_wmo_boxes; la_wmo_boxes_ctd(:,4)=0; 
la_wmo_boxes_argo=la_wmo_boxes; la_wmo_boxes_argo(:,2)=0; 
mkdir([owc_data_dir,filesep,'constants',filesep,'all']);
mkdir([owc_data_dir,filesep,'constants',filesep,'ctd']);
mkdir([owc_data_dir,filesep,'constants',filesep,'argo']);
				save([owc_data_dir,filesep,'constants',filesep,'all',filesep,'wmo_boxes'],'la_wmo_boxes');
la_wmo_boxes=la_wmo_boxes_ctd;	save([owc_data_dir,filesep,'constants',filesep,'ctd',filesep,'wmo_boxes'],'la_wmo_boxes');
la_wmo_boxes=la_wmo_boxes_argo;	save([owc_data_dir,filesep,'constants',filesep,'argo',filesep,'wmo_boxes'],'la_wmo_boxes');
% Then you can change the filenames to select, instead of fiddling inside matlab matrices

% ----- Check the reference data (i.e., make a lot of figures): -----
if true 
  for j=1:length(tartyp)	% Loop the three data types
    d=edir(tardir{j},'mat',0,0); %filer=cellstr(char(d.name));	% List of files
    hf=check_referencedata(d,true);
    close(hf);
  end % Loop the three data types
end 

