% A SCRIPT TO INGEST, QUICK-CHECK AND UPDATE LIST OF REFERENCE DATA
% FOR Argo DMQC.
%
% First, in a shell, download and unpack DMQC reference data:
% cd ~/Downloads/DMQC/
% lftp -u <user,password> ftp.ifremer.fr/coriolis/

clear all; close all;

init_dmqc;	% Get some necessary paths.

% Enter names of directories to ingest from:
refdir={'CTD_for_DMQC_2019V01_1','CTD_for_DMQC_2019V01_7','ARGO_for_DMQC_2019V03'};
refdir={'CTD_for_DMQC_2019V01_1','CTD_for_DMQC_2019V01_7','ARGO_for_DMQC_2020V01'};%May2020

% WMO-squares of chosen area to populate (fixed):
A=[1803 1804 1805	                         7802 7801 7800 1800 1801 1802 ...
   1703 1704 1705	7707 7706 7705           7702 7701 7700 1700 1701 1702 ...
   1603 1604 1605	     7606 7605 7604 7603 7602 7601 7600 1600 1601 1602 ...
			     7506 7505 7504 7503 7502 7501 7500 1500 1501 1502 ...
			     7406 7405 7404 7403 7402 7401 7400  ];

% Target directories (names of these subdirectories are given by OWC-toolbox):
tardir={[owc_data_dir,'climatology/historical_ctd'], ... 
	[owc_data_dir,'climatology/historical_bot'], ...
	[owc_data_dir,'climatology/argo_profiles']};
tartyp={'ctd_','bot_','argo_'};


% ----- Make an overview map of selected and ingested WMO-squares -----
figure(1); clf; set(gcf,'OuterPosition',get(0,'ScreenSize'));
wmolim=[]; for i=1:length(A), wmolim(i,:)=wmosquare(A(i)); end
mima(wmolim(:,3:4))+[-2 2];ans(2)=min(ans(2),89);
m_proj('Albers','lon',mima(wmolim(:,1:2))+[-2 2],'lat',ans);
m_grid; m_coast('color','k'); m_elev('contour','color',[.7 .7 .7]);
for i=1:length(A)
  [lo,la]=erect(wmolim(i,:),'t'); m_line(lo,la);
  sort([wmolim(i,1:2)+[-4 4],wmolim(i,1:2)]);lom=ans(2:3);
  sort([wmolim(i,3:4)+[-1 1],wmolim(i,3:4)]);lam=ans(2:3);
  htx=m_text(mean(lom),mean(lam),int2str(A(i))); 
  set(htx,'verticalalignment','middle','horizontalalignment','center','color','b');
end
% map will be printed after the content check below.



% ----- Transfer the reference data to the apropriate directories: -----
for j=1:length(refdir) % Loop reference-data directories
  if any(findstr('CTD',refdir{j}))
    for i=1:length(A)
      %['cp ',download_ref_data_dir,filesep,refdir{j},filesep,'ctd_',int2str(A(i)),'.mat ',tardir{1},filesep]
      msg=system(['cp ',download_ref_data_dir,filesep,refdir{j},filesep,'ctd_',int2str(A(i)),'.mat ',tardir{1},filesep]);
    end
  elseif any(findstr('ARGO',refdir{j}))
    for i=1:length(A)
      %['cp ',download_ref_data_dir,filesep,refdir{j},filesep,'argo_',int2str(A(i)),'.mat ',tardir{3},filesep]
      msg=system(['cp ',download_ref_data_dir,filesep,refdir{j},filesep,'argo_',int2str(A(i)),'.mat ',tardir{3},filesep]);
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
    if any(contains(filer,[tartyp{j},int2str(la_wmo_boxes(i,1)),'.mat']))
      la_wmo_boxes(i,1+j)=1;
      % Add text to the map for which data types exist in the WMOsq.
      wmosquare(la_wmo_boxes(i,1)); 
      la=linspace(ans(3),ans(4),10); la=la(1+j);
      htxex=m_text((ans(1)+ans(2))/2,la,tartyp{j}(1:end-1));
      set(htxex,'HorizontalAlignment','center');
    end
  end
end
print(gcf,'-depsc',[owc_data_dir,filesep,'constants',filesep,'selected_wmos_map.eps']);
save([owc_data_dir,filesep,'constants',filesep,'wmo_boxes'],'la_wmo_boxes');


% ----- Check the reference data (make a lot of figures): -----
if logical(1) 
  %figure(2);set(gcf,'OuterPosition',[1 385 1026 960]);
  %figure(2);set(gcf,'OuterPosition',[1 385 1426 960]);
  figure(2);set(gcf,'OuterPosition',get(0,'ScreenSize'));
  for j=1:length(tartyp)	% Loop the three data types
    d=edir(tardir{j},'mat',0,0); %filer=cellstr(char(d.name));	% List of files
    for i=1:length(d)	% Loop all mat-files
      load(d{i});
      long>180;long(ans)=long(ans)-360;
      % Check positions:
      wmosq=str2num(d{i}(end-7:end-4));
      if unique(findwmo(long,lat)) ~= wmosq % Check positions
	warning(['There are ',upper(tartyp{j}(1:end-1)),'-positions outside the WMO-square ',int2str(wmosq),' !'])
      end
      % Check for zero sal and temp _and_ sal:
      sal==0 | sal==0 & temp==0;
      if any(ans,'all')
	warning(['There are zero S&T, or just S, in ',upper(tartyp{j}(1:end-1)),...
		 ' data WMO-square ',int2str(wmosq),' ! (removed)'])
	sal(ans)=NaN; temp(ans)=NaN; % Useless in both cases
	save(d{i},'sal','temp','-append');
      end      
      % Check for incomplete sal and temp pairs:
      isnan(sal) & ~isnan(temp) | ~isnan(sal) & isnan(temp);
      if any(ans,'all')
	warning(['There are incomplete S&T pairs in ',upper(tartyp{j}(1:end-1)),...
		 ' data WMO-square ',int2str(wmosq),' !'])
	sal(ans)=NaN; temp(ans)=NaN;
      end      
      
      % Figure:
      figure(2);clf;

      % Map with colors for time:
      subplot 231; 
      %title([upper(tartyp{j}(1:end-1)),'-reference data in WMO-square ',int2str(wmosq)]);
      title(['Time span ',datestr(datenum(min(dates)/1e10,0,0),12),'-',datestr(datenum(max(dates)/1e10,0,0),12)]);
      lim=wmosquare(wmosq);
      m_proj('Albers','lon',lim(1:2)+[-4 4],'lat',[lim(3)-1 min([90 lim(4)+1])]);
      m_grid; m_coast('color','k'); m_elev('contour','color',[.7 .7 .7]);
      [lo,la]=erect(lim,'t'); hw=m_line(lo,la,'color','b');
      [ans,IA]=sort(dates);
      %hp=m_line(long,lat,'linestyle','none','marker','.','color','b');
      hp=m_scatter(long(IA),lat(IA),20,jet(size(pres,2)));set(hp,'marker','.');
      hl=legend([hw hp(1)],'WMO-square','Reference data (by time)','location','northwestoutside');

      % Map with colors for latitude (as for the rest):
      subplot 232; 
      title([upper(tartyp{j}(1:end-1)),'-reference data in WMO-square ',int2str(wmosq)]);
      lim=wmosquare(wmosq);
      %m_proj('Albers','lon',lim(1:2)+[-4 4],'lat',[lim(3)-1 min([90 lim(4)+1])]);
      m_grid; m_coast('color','k'); m_elev('contour','color',[.7 .7 .7]);
      [lo,la]=erect(lim,'t'); hw=m_line(lo,la,'color','b');
      [ans,IA]=sort(lat);
      %hp=m_line(long,lat,'linestyle','none','marker','.','color','b');
      hp=m_scatter(long(IA),lat(IA),20,jet(size(pres,2)));set(hp,'marker','.');
      hl=legend([hw hp(1)],'WMO-square','Reference data (by latitude)','location','northeastoutside');

      % TS-diagram:
      subplot 234
      tsdiagrm(mima(sal),mima(temp),0);
      title('Reference data (by latitude)');
      hTS=line(sal(:,IA),temp(:,IA),'marker','.','linestyle','none');
      set(hTS,{'color'},num2cell(jet(size(pres,2)),2));
      %set(gca,'xlim',[31.5 35.5],'ylim',[-2 16]); % Lock for specific presentation
      set(hTS,'clipping','off');
      
      % Profiles:
      subplot 235
      hT=line(pres(:,IA),temp(:,IA)); 
      view([90 90]);grid;set(gca,'yaxislocation','right');
      set(hT,{'color'},num2cell(jet(size(pres,2)),2));
      xlabel Pressure; ylabel Temperature
      %set(gca,'xlim',[0 3500],'ylim',[-2 16]);  % Lock for specific presentation
      set(hT,'clipping','off');
      subplot 236
      hS=plot(pres(:,IA),sal(:,IA)); 
      view([90 90]);grid;set(gca,'yaxislocation','right');
      set(hS,{'color'},num2cell(jet(size(pres,2)),2));
      xlabel Pressure; ylabel Salinity
      %set(gca,'xlim',[0 3500],'ylim',[31.5 35.5]);  % Lock for specific presentation
      set(hS,'clipping','off');
      %
      print(gcf,'-depsc',[d{i}(1:end-4),'.eps']); 
    end
  end
end


%set([hTS,hT,hS],'',
      
