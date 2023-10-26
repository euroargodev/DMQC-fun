% PLOT_PROFILES makes overview plots of the positions and T&S data.
%
% PREPARE_FLOATS uses this to make the five panel plots with map,
% TS-diagram and profile plots of T, S and density: the raw input; with
% flagging by the automated tests prior to visual control marked; and
% finally with the cleaned up data only.
%
% DMQC-fun v0.9.
% J. Even Ø. Nilsen, Ingrid M. Angel-Benavides, Birgit Klein, Malgorzata Merchel, and Kjell Arne Mork.
% Last updated: Tue Oct 24 22:47:55 2023 by jan.even.oeie.nilsen@hi.no


%LONG(LONG>180)=LONG(LONG>180)-360;
% 28.05.2020: Corrected the erroneous change of LONG to negative values west of zero meridian:
lon=LONG; lon(lon>180)=lon(lon>180)-360;

figure(str2num(char(float_names{I})));clf;
%set(gcf,'OuterPosition',get(0,'ScreenSize'));
shape=[1 1 1026 960];
set(gcf,'units','points','innerposition',shape,'paperunits','points','paperposition',shape,'PaperSize',shape(3:4),'PaperPositionMode','manual','RendererMode','manual','Renderer','opengl');

[hT,hS,hD]=deal([]);

% Map:
a_map=subplot(2,3,1); 
%title([upper(tartyp{j}(1:end-1)),'-reference data in WMO-square ',int2str(wmosq)]);
if exist('nbt'), title([nbt,' at float number ',float_names{I}]); 
else,            title(float_names{I});					end
lim=[mima(lon) mima(LAT)];
m_proj('Albers','lon',lim(1:2)+[-4 4],'lat',lim(3:4)+[-1 1])
m_elev('contour',[-5000:250:0],'color',[.7 .7 .7]);
m_coast('color','k'); 
m_grid; 
wmos=unique(findwmo(lon,LAT));wmos=wmos(~isnan(wmos));
for ii=1:length(wmos)
  wmolim=wmosquare(wmos(ii)); [lo,la]=erect(wmolim,'t'); m_line(lo,la);
  sort([lim(1:2)+[-4 4],wmolim(1:2)]);lom=ans(2:3);
  sort([lim(3:4)+[-1 1],wmolim(3:4)]);lam=ans(2:3);
  htx=m_text(mean(lom),mean(lam),int2str(wmos(ii)));
  set(htx,'verticalalignment','middle','horizontalalignment','center','color','b');
end
%[lo,la]=erect(lim,'t'); hw=m_line(lo,la,'color','b');
IA=1:length(lon); %[ans,IA]=sort(LAT);
hp=m_scatter(lon,LAT,5,jet(n));
set(hp,'marker','*');
%ht=m_text(lon(10:10:end),LAT(10:10:end),int2str(PROFILE_NO(10:10:end)'));
[1 10:10:n n]; ht=m_text(lon(ans),LAT(ans),int2str(CYCLE_NUMBER(ans)'));
if exist('jnb'), hpnb=m_line(lon(jnb),LAT(jnb)); set(hpnb,'marker','o','color','k','linestyle','none'); end
%hl=legend([hw hp(1)],'WMO-square','Reference data','location','northwestoutside');
hcb=colorbar('southoutside');colormap(jet(n));caxis(mima(CYCLE_NUMBER)+[-.01 .01]);%caxis([0 n]);
%xlabel(hcb,['Cycle number (',datestr(datenum(DATES(1),0,0),12),'-',datestr(datenum(DATES(end),0,0),12),')']);
xlabel(hcb,['Cycle number (',datestr(time(1),12),'-',datestr(time(end),12),')']);
% Sea ice concentration if lost position:
intpos=POSqco=='8'; 
if any(intpos)
  groups(intpos); ans(~intpos)=0; [~,IA]=unique(ans,'stable');
  for is=2:length(IA) % loop periods of missing positions
    esic(lon(IA(is)-1),LAT(IA(is)-1),time(IA(is)),2);
  end % Plot SIC on map near last valid position at day of disappearance
end

% TS-diagram:
a_TS=subplot(2,3,2);
if ~all(isnan(PSAL_ADJUSTED),'all') & ~all(isnan(TEMP_ADJUSTED),'all')
  tsdiagrm(mima(PSAL_ADJUSTED),mima(TEMP_ADJUSTED),0);
  hTS=line(PSAL_ADJUSTED,TEMP_ADJUSTED,'marker','.','linestyle','none');
  if exist('snb'), hTSnb=line(PSAL_ADJUSTED(snb),TEMP_ADJUSTED(snb),'color','k','linewidth',2,'linestyle','none','marker','o'); end
  if exist('tnb'), hTSnb=line(PSAL_ADJUSTED(tnb),TEMP_ADJUSTED(tnb),'color','k','linewidth',2,'linestyle','none','marker','o'); end
  if exist('inb'), hTSnb=line(PSAL_ADJUSTED(inb),TEMP_ADJUSTED(inb),'color','k','linewidth',2,'linestyle','none','marker','o'); end
  set(hTS,{'color'},num2cell(jet(size(PRES_ADJUSTED,2)),2));
else
  text(mean(xlim),mean(ylim),'No T/S data to show','HorizontalAlignment','center'); axis off;
end

% Profiles:

aT=subplot(2,3,4);
if ~all(isnan(TEMP_ADJUSTED),'all')
  hT=line(PRES_ADJUSTED,TEMP_ADJUSTED); 
  if exist('tnb'), hTnb=line(PRES_ADJUSTED(tnb),TEMP_ADJUSTED(tnb),'color','k','linewidth',2,'linestyle','none','marker','o'); end
  if exist('inb'), hTnb=line(PRES_ADJUSTED(inb),TEMP_ADJUSTED(inb),'color','k','linewidth',2,'linestyle','none','marker','o'); end
  view([90 90]);grid;set(gca,'yaxislocation','right','box','on');
  set(hT,{'color'},num2cell(jet(size(PRES_ADJUSTED,2)),2));
  xlabel Pressure; ylabel Temperature
else
  text(mean(xlim),mean(ylim),'No T data to show','HorizontalAlignment','center'); axis off; aT=[];
end

aS=subplot(2,3,5);
if ~all(isnan(PSAL_ADJUSTED),'all')
  hS=line(PRES_ADJUSTED,PSAL_ADJUSTED); 
  if exist('snb'), hSnb=line(PRES_ADJUSTED(snb),PSAL_ADJUSTED(snb),'color','k','linewidth',2,'linestyle','none','marker','o'); end
  if exist('inb'), hSnb=line(PRES_ADJUSTED(inb),PSAL_ADJUSTED(inb),'color','k','linewidth',2,'linestyle','none','marker','o'); end
  view([90 90]);grid;set(gca,'yaxislocation','right','box','on');
  set(hS,{'color'},num2cell(jet(size(PRES_ADJUSTED,2)),2));
  xlabel Pressure; ylabel Salinity
else
  text(mean(xlim),mean(ylim),'No S data to show','HorizontalAlignment','center'); axis off; aS=[];
end

aD=subplot(2,3,6);
if ~all(isnan(DENS),'all')
  hD=line(PRES_ADJUSTED,DENS); 
  if exist('snb'), hDnb=line(PRES_ADJUSTED(snb),DENS(snb),'color','k','linewidth',2,'linestyle','none','marker','o'); end
  if exist('tnb'), hDnb=line(PRES_ADJUSTED(tnb),DENS(tnb),'color','k','linewidth',2,'linestyle','none','marker','o'); end
  if exist('pnb'), hDnb=line(PRES_ADJUSTED(pnb),DENS(pnb),'color','k','linewidth',2,'linestyle','none','marker','o'); end
  if exist('inb'), hDnb=line(PRES_ADJUSTED(inb),DENS(inb),'color','k','linewidth',2,'linestyle','none','marker','o'); end
  view([90 90]);grid;set(gca,'yaxislocation','right','box','on');
  set(hD,{'color'},num2cell(jet(size(PRES_ADJUSTED,2)),2));
  xlabel Pressure; ylabel('Potential Density')
else
  text(mean(xlim),mean(ylim),'No density data to show','HorizontalAlignment','center'); axis off; aD=[];
end

% Reposition the two upper panels:
get(a_map,'position'); set(a_map,'position',ans+[0 -.07 .15 .12]);
get(a_TS,'position'); set(a_TS,'position',ans+[.18 -.05 .1 .1]); 

% Set same pressure limits on all:
mima(get([hT;hS;hD],'xdata')); set([aT,aS,aD],'xlim',ans+[-10 10]);


