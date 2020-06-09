% plotting for load_floats.m

%LONG(LONG>180)=LONG(LONG>180)-360;
% 28.05.2020: Corrected the erroneous change of LONG to negative values west of zero meridian:
lon=LONG; lon(lon>180)=lon(lon>180)-360;

figure(str2num(char(float_names{I})));clf;
%set(gcf,'OuterPosition',get(0,'ScreenSize'));
set(gcf,'OuterPosition',[1 385 1026 960]);

% Map:
subplot 231; 
%title([upper(tartyp{j}(1:end-1)),'-reference data in WMO-square ',int2str(wmosq)]);
if exist('nbt'), title([nbt,' at float number ',float_names{I}]); end
lim=[mima(lon) mima(LAT)];
m_proj('Albers','lon',lim(1:2)+[-4 4],'lat',lim(3:4)+[-1 1])
m_grid; m_coast('color','k'); m_elev('contour','color',[.7 .7 .7]);
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
hp=m_scatter(lon,LAT,5,jet(n));set(hp,'marker','*');
ht=m_text(lon(10:10:end),LAT(10:10:end),int2str(PROFILE_NO(10:10:end)'));
if exist('jnb'), hpnb=m_line(lon(jnb),LAT(jnb)); set(hpnb,'marker','o','color','k','linestyle','none'); end
%hl=legend([hw hp(1)],'WMO-square','Reference data','location','northwestoutside');
hcb=colorbar('southoutside');colormap(jet(n));caxis([0 n]);
xlabel(hcb,['Profile number (',datestr(datenum(DATES(1),0,0),12),'-',datestr(datenum(DATES(end),0,0),12),')']);

% TS-diagram:
subplot 232
tsdiagrm(mima(SAL),mima(TEMP),0);
hTS=line(SAL,TEMP,'marker','.','linestyle','none');
if exist('snb'), hTSnb=line(SAL(snb),TEMP(snb),'color','k','linewidth',2,'linestyle','none','marker','o'); end
if exist('tnb'), hTSnb=line(SAL(tnb),TEMP(tnb),'color','k','linewidth',2,'linestyle','none','marker','o'); end
if exist('inb'), hTSnb=line(SAL(inb),TEMP(inb),'color','k','linewidth',2,'linestyle','none','marker','o'); end
set(hTS,{'color'},num2cell(jet(size(PRES,2)),2));

% Profiles:
subplot 234
hT=line(PRES,TEMP); 
if exist('tnb'), hTnb=line(PRES(tnb),TEMP(tnb),'color','k','linewidth',2,'linestyle','none','marker','o'); end
if exist('inb'), hTnb=line(PRES(inb),TEMP(inb),'color','k','linewidth',2,'linestyle','none','marker','o'); end
view([90 90]);grid;set(gca,'yaxislocation','right');
set(hT,{'color'},num2cell(jet(size(PRES,2)),2));
xlabel Pressure; ylabel Temperature

subplot 235
hS=line(PRES,SAL); 
if exist('snb'), hSnb=line(PRES(snb),SAL(snb),'color','k','linewidth',2,'linestyle','none','marker','o'); end
if exist('inb'), hSnb=line(PRES(inb),SAL(inb),'color','k','linewidth',2,'linestyle','none','marker','o'); end
view([90 90]);grid;set(gca,'yaxislocation','right');
set(hS,{'color'},num2cell(jet(size(PRES,2)),2));
xlabel Pressure; ylabel Salinity

subplot 236
hD=line(PRES,DENS); 
if exist('snb'), hDnb=line(PRES(snb),DENS(snb),'color','k','linewidth',2,'linestyle','none','marker','o'); end
if exist('tnb'), hDnb=line(PRES(tnb),DENS(tnb),'color','k','linewidth',2,'linestyle','none','marker','o'); end
if exist('pnb'), hDnb=line(PRES(pnb),DENS(pnb),'color','k','linewidth',2,'linestyle','none','marker','o'); end
if exist('inb'), hDnb=line(PRES(inb),DENS(inb),'color','k','linewidth',2,'linestyle','none','marker','o'); end
view([90 90]);grid;set(gca,'yaxislocation','right');
set(hD,{'color'},num2cell(jet(size(PRES,2)),2));
xlabel Pressure; ylabel Density

%try 
%  set([hTSnb,hTnb,hSnb,hDnb],'color','k','linewidth',2,'linestyle','none','marker','o')
%end

%delete([hTSnb,hTnb,hSnb,hpnb]); % For a clean plot


% Have to change more here: 
% Better marking of inversions/non-monotonic parts;
% Plot density too. (PDENS?)
