function hf = check_referencedata(file,wrte)
% CHECK_REFERENCEDATA	Plots and checks Argo RDB files
% 
% hf = check_referencedata(file)
%
% file  = cellstring list of filenames or a single character file
%	  name. To enter a list of all mat files in a directory you
%	  can use EDIR.  
% wrte  = logical whether to make corrections in the input files
%	  (false by default).
% 
% hf    = handle to figure window.
% 
% Other outputs are png-files of figures alongside the input
% files. Furthermore, the file content may be altered (if wrte): In the
% case there are no S&T, or just S, data these points are removed as
% they are useless for OWC.  Other discrepancises are warned in command
% line output.
%
% See also EDIR FINDWMO WMOSQUARE

% Last updated: Thu Oct 26 15:27:54 2023 by jan.even.oeie.nilsen@hi.no

error(nargchk(1,2,nargin));
if nargin<2 | isempty(wrte),   wrte=false;   end

file=cellstr(file);

% Prepare figure window:
hf=figure; set(hf,'OuterPosition',get(0,'ScreenSize'));
shape=get(hf,'innerposition');
set(hf,'units','points','innerposition',shape,'paperunits','points','paperposition',shape,...
       'PaperSize',shape(3:4),'PaperPositionMode','manual','RendererMode','manual','Renderer','opengl');


for i=1:length(file)		% Loop all mat-files
  
  if endsWith(file{i},'.mat'), file{i}=replace(file{i},'.mat',''); end
  %split(file{i},{'_'});
  split(file{i},{filesep,'_'});  
  tartyp=ans{end-1};
  wmosq=str2num(ans{end});
  load(file{i});
  D=size(pres);
  long>180;long(ans)=long(ans)-360;

  % Check positions:
  findwmo(long,lat)~=wmosq;
  if any(ans,'all') % Check positions
    J=unique(find(ans));
    disp(['Profiles ',zipnumstr(J),...
	  ' in ',upper(tartyp),...
	  ' data for ',int2str(wmosq),...
	  ' are outside the WMO-square!']);
  end

  % Check for zero sal and temp _and_ sal:
  sal==0 | sal==0 & temp==0;
  if any(ans,'all')
    [I,J] = ind2sub(D,find(ans)); I=unique(I); J=unique(J);
    msg=['There are zero S&T, or just S, in ',upper(tartyp),...
	 ' data in WMO-square ',int2str(wmosq),...
	 ' typically on rows ',zipnumstr(I),...
	 ' in profiles ',zipnumstr(J),...
	 ' from source ',snippet(unique(source(J)),';'),...
	 ' ! (removed)'];
    if exist('qclevel'), msg=[msg,' The qclevels are ',snippet(unique(qclevel(J)),';'),'.']; end
    disp(msg)
    sal(ans)=NaN; temp(ans)=NaN; % Useless in both cases
    if wrte, save(file{i},'sal','temp','-append'); end
  end      

  % Check for incomplete sal and temp pairs:
  isnan(sal) & ~isnan(temp) | ~isnan(sal) & isnan(temp);
  if any(ans,'all')
    [I,J] = ind2sub(D,find(ans)); I=unique(I); J=unique(J);
    msg=['There are incomplete S&T pairs in ',upper(tartyp),...
	 ' data in WMO-square ',int2str(wmosq),...
	 ' typically on rows ',zipnumstr(I),...
	 ' in profiles ',zipnumstr(J),...
	 ' from source ',snippet(unique(source(J)),';'),...
	 ' !'];
    if exist('qclevel'), msg=[msg,' The qclevels are ',snippet(unique(qclevel(J)),';'),'.']; end
    disp(msg)
    sal(ans)=NaN; temp(ans)=NaN;
    if wrte, save(file{i},'sal','temp','-append'); end
  end      
  
  % Goto and clear figure:
  figure(hf); clf;
  
  % Map with colors for time:
  a_map1=subplot('position',[0.1 0.5 0.35 0.4]); 
  title(['Time span ',datestr(rdtime(min(dates)),12),'-',datestr(rdtime(max(dates)),12)]);
  lim=wmosquare(wmosq);
  m_proj('Albers','lon',lim(1:2)+[-4 4],'lat',[lim(3)-1 min([90 lim(4)+1])]);
  m_grid; m_coast('color','k'); m_elev('contour','color',[.7 .7 .7]);
  [lo,la]=erect(lim,'t'); hw=m_line(lo,la,'color','b');
  [ans,IA]=sort(rdtime(dates));
  hp=m_scatter(long(IA),lat(IA),20,jet(size(pres,2)));set(hp,'marker','.');
  hl=legend([hw hp(1)],'WMO-square','Reference data (by time)','location','northwestoutside');
  
  % Map with colors for latitude (as for the rest):
  a_map2=subplot('position',[0.55 0.5 0.35 0.4]); 
  title([upper(tartyp),'-reference data in WMO-square ',int2str(wmosq)]);
  lim=wmosquare(wmosq);
  m_grid; m_coast('color','k'); m_elev('contour','color',[.7 .7 .7]);
  [lo,la]=erect(lim,'t'); hw=m_line(lo,la,'color','b');
  [ans,IA]=sort(lat);
  hp=m_scatter(long(IA),lat(IA),20,jet(size(pres,2)));set(hp,'marker','.');
  hl=legend([hw hp(1)],'WMO-square','Reference data (by latitude)','location','northeastoutside');
  
  % TS-diagram:
  subplot 234
  tsdiagrm(mima(sal),mima(temp),0);
  title('Reference data (by latitude)');
  hTS=line(sal(:,IA),temp(:,IA),'marker','.','linestyle','none');
  set(hTS,{'color'},num2cell(jet(size(pres,2)),2));
  set(hTS,'clipping','off');
      
  % Profiles:
  subplot 235
  hT=line(pres(:,IA),temp(:,IA)); 
  view([90 90]);grid;set(gca,'yaxislocation','right');
  set(hT,{'color'},num2cell(jet(size(pres,2)),2));
  xlabel Pressure; ylabel Temperature
  set(hT,'clipping','off');
  subplot 236
  hS=plot(pres(:,IA),sal(:,IA)); 
  view([90 90]);grid;set(gca,'yaxislocation','right');
  set(hS,{'color'},num2cell(jet(size(pres,2)),2));
  xlabel Pressure; ylabel Salinity
  set(hS,'clipping','off');
      
  % Print to png-file alongside the mat-file: 
  print(hf,'-depsc',[file{i},'.eps']); 

  clear dates lat long pres ptmp qclevel sal source temp
      
end % Loop all mat-files
