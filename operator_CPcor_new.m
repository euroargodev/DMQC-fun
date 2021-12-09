% OPERATOR_CPCOR_NEW helps find the best CPcor_new
% by J. Even Ø. Nilsen, Ingrid Angel, Birgit Klein, and Kjell Arne Mork.
% with contribution from Cécile Cabanes.
% DMQC-fun v0.9.3, jan.even.oeie.nilsen@hi.no.
%
% Not a beautiful solution, but it works. However, it can be
% considered optional since it is not entirely necessary (and only
% relevant for deep floats), and really not user friendly as of now.
%
% You can use this script to find operator's CPcor_new and save it to
% your working directory for your float. You only need to do this
% once, so you can use any apporach, as long as you save the results
% as done at the end of this script. PREPARE_FLOATS will then use the
% saved results every time you do DMQC. 
% 
% This script is bigger than it needs to be, but this is put together to
% explore how to best find a better CPcor-new. With the settings as is
% and the decision 'king', you will get the CPcor_new according to
% Cabanes/King.
%
% The code for this default method is from the toolbox DM_CPcor by
% Cécile Cabanes (https://github.com/ArgoDMQC/DM_CPcor/releases). It has
% been adapted here (courtesy of Cécile Cabanes) for simpler integration
% with the DMQC-fun system.  You can of course use Cécile Cabanes'
% toolbox instead of this script, and then save the CPcor_new, ctd_data,
% and other parameters for PREPARE_FLOATS to use, in the same manner as
% is done in the save section at the end of this script.
%
% Initial original code is from argo_test_dypArgo, by Kjell Arne Mork,
% incorporated into DMQC-fun and furter developed by J.Even.Ø.Nilsen.
%
% Needs CTD data put in a 'CTD' directory of the float's DMQC directory.
% Results and figures are put in 'CPcor' directory of the float's DMQC directory.
% See WORK_LOG for explanation of when to use.
%
% OPTIONAL TOOLBOXES:
%	seawater
%	ctdbase (gctd.m) by kjell.arne.mork@hi.no can be used to read
%	seabird cnv (or kva) files while running this script. 
%
% MANUAL INPUT OF CTD DATA: Put your cruise CTD data in a struct called
% 'ctd' with fields
%
%               t: [1xN double]
%             lat: [1xN double]
%             lon: [1xN double]
%         station: [1xN double]
%            temp: [MxN double]
%            salt: [MxN double]
%               z: [Mx1 double]
%	       
% where t is in Matlab serial date number (see DATENUM); station is the
% station number from the cruise; and z is pressure/depth.  You will
% enter which station (column) to use in the table below. Make a
% subdirectory called 'CTD' in <my_working_dir>/DMQC/<float_name> and
% save this object in a file called 'CTD_<float_name>.mat'.
%
% THE DECISION-MAKING PROCESS: This script is run after PREPARE_FLOATS,
% as it uses early stage DMQC'd data and some other input from
% there. You will have to maintain the list of float names and
% corresponding cruise and CTD station numbers in this script,
% below. When you have found the station to use and decided on method of
% CPcor_new determination, you run for the last time, and re-run
% PREPARE_FLOATS.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; 
init_dmqc; % Paths and filenames. Set one float at a time when
           % considering new ctd data and decisions.

scr=get(0,'screensize'); figz=ceil(scr*.6);
CPcor_SBE = -9.57e-8;
tool='sw'; % 'sw' or 'gsw' % GSW DOES NOT WORK IN THIS SCRIPT! DO NOT USE!
cecile=true; % If cecile's results are to be added

for I=1:length(float_names)	% Loop floats
  disp(['Finding operator''s CPcor_new for float ',float_names{I},'.']);
  close all
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EDIT THIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % THE LIST OF MATCHING FILES AND DECISIONS YOU MAKE ('rec', 'var', 'sqr', 'king', or 'gain'):
  % np is the Argo profile number (the first in some cases is too short).  
  switch float_names{I}
   case '6903556', ctd_cruise='2020608', ctd_file='4-2020-1172-8_tkt_kva.txt';	ctd_station='435';	np=1; decision='king';%√!var=par
   case '6903557', ctd_cruise='2019205', ctd_file='4-2019-1019-5_tkt_kva.txt';	ctd_station='501';	np=1; decision='king';%√sqr%502
   case '6903558', ctd_cruise='2019205', ctd_file='4-2019-1019-5_tkt_kva.txt';	ctd_station='547';	np=1; decision='king';%√sqr/var?%548
   case '6903571', ctd_cruise='2020608', ctd_file='4-2020-1172-8_tkt_kva.txt';	ctd_station='425';	np=1; decision='king';%√sqr/var?
   case '6903572', ctd_cruise='2020616', ctd_file='4-2020-1172-16_tkt_kva.txt';	ctd_station='684';	np=1; decision='king';%√!var=par
   case '6903573', ctd_cruise='2020616', ctd_file='4-2020-1172-16_tkt_kva.txt';	ctd_station='684';	np=1; decision='king';%714grunn%684%√!var
   case '6903580', ctd_cruise='2021205', ctd_file='T21JH205.KVA';		ctd_station='358';	np=2; decision='king';%√sqr
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % ------- Subdirectories of float DMQC results (automatic): ----------
  cpcordir = [my_working_dir,'DMQC',filesep,float_names{I},filesep,'CPcor']; if ~exist(cpcordir,'dir'), mkdir(cpcordir); end
  ctddir   = [my_working_dir,'DMQC',filesep,float_names{I},filesep,'CTD']; if ~exist(ctddir,'dir'), mkdir(ctddir); end
  
  % ------- CTD data: --------------------------------------------------
  ctdfile=[ctddir,filesep,'CTD_',float_names{I},'.mat']; 
  if exist(ctdfile,'file'), load(ctdfile);	% If already done, or constructed in other way, just load it. 
  else, ctd=gctd; save(ctdfile,'ctd');		% If not, read reference data from files on this list, using KAM's gctd GUI.
  end
  % which column to use in ctd matrix:
  ind=find(ctd.station==str2num(ctd_station)); 

  % ------- Argo float data: -------------------------------------------
  % Read from the float data prepared for OWC by PREPARE_FLOATS (from Rfiles
  % and with initial DMQC but no OWC done on it: ). First profile is first
  load(outfiles{I});
  flt.pres=PRES;
  flt.temp=TEMP;
  flt.psal=PSAL; % The one just saved before correcting CPcor.
  flt.ptmp=PTMP;
  flt.lon=LONG;
  flt.lat=LAT;
  flt.juld=datenum(DATES,1,1); 
  % [√] Has flt.juld taken into account the REFERENCE_DATE? Yes. 
  % [√] Why 2? KAM: det kan være at første profil ikke var så dyp (det kan
  % være at man velger å ta første profil kort tid etter utsetting slik at
  % bøyen ikke rekker å gå helt ned). Jeg mener at generelt bør første
  % profil velges hvis det er en dyp profil.
  argo_p=flt.pres(:,np);
  argo_t=flt.temp(:,np);
  argo_s=flt.psal(:,np);
  pt=flt.ptmp(:,np); %pt=sw_ptmp(s,t,p,0); % made by same method in PREPARE_FLOATS

  % ------- Positions: -------------------------------------------------
  % Common position for GSW to use, for consistency:
  lon=mean([flt.lon(np),ctd.lon(ind)]); lat=mean([flt.lat(np),ctd.lat(ind)]);
  lon<0; lon(ans)=lon(ans)+360;
  % Strings about time and place for Argo and CTD:
  flt.lon>180; flt.lon(ans)=flt.lon(ans)-360; ctd.lon>180; ctd.lon(ans)=ctd.lon(ans)-360;
  % floatstr=['Argo: ' datestr(flt.juld(np)) '; ' num2str(flt.lon(np),'%10.3f') ' , ' num2str(flt.lat(np),'%10.3f')]
  % leg.ctd=['Ship: ' datestr(ctd.t(ind)) '; ' num2str(ctd.lon(ind),'%10.3f') ' , ' num2str(ctd.lat(ind),'%10.3f')]
  leg.float=[datestr(flt.juld(np)),', ',num2str(flt.lon(np),'%10.3f'),'\circE, ',num2str(flt.lat(np),'%10.3f'),'\circN, Argo ',float_names{I},', cycle ',int2str(np)]
  leg.ctd=[datestr(ctd.t(ind)),', ',num2str(ctd.lon(ind),'%10.3f'),'\circE, ',num2str(ctd.lat(ind),'%10.3f'),'\circN, Cruise ',ctd_cruise,', station ',ctd_station]

  % ------- CTD-profiles to use: ---------------------------------------
  ctd_z    = ctd.z; % for PREPARE_FLOATS
  ctd_salt = ctd.salt(:,ind);
  ctd_temp = ctd.temp(:,ind);
  switch tool
   case 'sw'
    ctd_ptmp   = sw_ptmp(ctd_salt,ctd_temp,ctd_z,0);
   case 'gsw'
    ctd.SA     = gsw_SA_from_SP(ctd_salt,ctd_z,lon,lat); % (SP,p,long,lat)
    ctd_ptmp   = gsw_pt_from_t(ctd.SA,ctd_temp,ctd_z,0); % (SA,t,p,p_ref)
  end

  % ------- Finding CPcor_new: -----------------------------------------
  % The steps to re-compute salinity with CPcor_new are as follows.
  % (a). Fill PRES_ADJUSTED and TEMP_ADJUSTED [] ??????  Follow the same
  % 'D' mode procedures as for the 2000-dbar floats, as in Sections 3.3
  % and 3.4.
  % (b). Compute original conductivity, Co.  This is done by using PRES,
  % TEMP, PSAL. For example, if using the Gibbs-SeaWater (GSW)
  % Oceanographic Matlab Toolbox, then Co = gsw_C_from_SP (PSAL, TEMP, PRES).
  Co=NaN(size(argo_p));
  ii=find(~isnan(argo_p)&~isnan(argo_s)&~isnan(argo_t));	% All parameters non-NaN because gsw_ can't take any NaNs! :-( 
  Co(ii)=gsw_C_from_SP(argo_s(ii),argo_t(ii),argo_p(ii)); %(SP,t,p)	%%%% cond_argo
  d = 3.25e-6;						%%%% CTcor_SBE
  % CPcor_SBE = -9.57e-08 dbar-1,
  % CPcor_new = -12.5e-08 dbar-1 for SBE-61 data,
  % CPcor_new = -13.5e-08 dbar-1 for Deep SBE-41CP data, or
  CPcor_SBE = -9.57e-8;
  % Calculate psal and ptemp profile for many CPcor_new:
  clear CPcor_new s_new Cnew pt_new ds_var ds_sqr ds_glatt winsize
  Ncp=150;
  for i=1:Ncp
    CPcor_new.f(i) = (-20+i/10)*1.0e-8;
    switch tool
     case 'sw'
      s_new.f(:,i) = argo_sbe_CPcor(argo_p,argo_t,argo_s,CPcor_new.f(i));
      pt_new.f(:,i) = sw_ptmp(s_new.f(:,i),argo_t,argo_p,0);
     case 'gsw'
      % (c). Compute new conductivity, Cnew.
      % Cnew = Co*(1 + d*TEMP + CPcor_SBE*PRES) / (1 + d*TEMP_ADJUSTED + CPcor_new*PRES_ADJUSTED);
      % where d = 3.25e-06,
      Cnew.f(:,i) = Co .* (1 + d*argo_t + CPcor_SBE*argo_p) ./ (1 + d*argo_p + CPcor_new.f(i)*argo_p); %%%%cond_argo_DM=cond_ZERO/b1=cond_argo*a1/b1
      % (d). Compute new salinity, PSAL_ADJUSTED_Cnew.  This is done by
      % using PRES_ADJUSTED, TEMP_ADJUSTED, and Cnew. For example, if
      % using the Gibbs-SeaWater (GSW) Oceanographic Matlab Toolbox, then
      % PSAL_ADJUSTED_Cnew = gsw_SP_from_C (Cnew, TEMP_ADJUSTED, PRES_ADJUSTED).
      s_new.f(:,i) = gsw_SP_from_C(Cnew.f(:,i),argo_t,argo_p); % (C,t,p)			%%%%psal_argo_corr
      flt.SA= gsw_SA_from_SP(s_new.f(:,i),argo_p,lon,lat);  % (SP,p,long,lat)	%%%%
      pt_new.f(:,i) = gsw_pt_from_t(flt.SA,argo_t,argo_p,0);  % (SA,t,p,p_ref)
    end
    % Calculate new profiles for this i-th CPcor_new:
    indt=find(pt_new.f(:,i)<-0.5); %indt=find(argo_p>1500);	% limit to cold/deep waters
    plim(i)=min(argo_p(indt)); % The depth which we consider. 
    inda=find(~isnan(ctd_ptmp) & ~isnan(ctd_salt));	% non-nan values only
    ctdsalt_a = interp1(ctd_ptmp(inda),ctd_salt(inda),pt_new.f(indt,i));
    %ctdsalt_a = interp1(ctd_z(inda),ctd_salt(inda),argo_p(indt));
    ds = s_new.f(indt,i) - ctdsalt_a;				% the deviations in salt
    %[ds_glatt, winsize] = smoothdata(ds,'movmedian',15,'omitnan');
    [ds_glatt, winsize] = smoothdata(ds,'movmean',15,'omitnan');% smooth deviations
    ds_var(i)=nanvar(ds_glatt);					% Variance of deviations in profile
    ds_sqr(i)=nansum(ds_glatt.*ds_glatt);			% Sum of squared deviations in profile
  end
  % Find the smallest deviation (two methods):
  [dsv_m,indvar]=nanmin(ds_var);
  [dsv_q,indsqr]=nanmin(ds_sqr);

  % ------- Find the respective CPcor: --------------------------------- 
  CPcor_old=CPcor_SBE;
  switch ctdmodel % loaded with float data
   case 'SBE41CP', CPcor_new.rec = -13.5e-8;
   case 'SBE61',   CPcor_new.rec = -12.5e-8;
  end
  CPcor_new.var=CPcor_new.f(indvar); % ds1=nanmean(s_new.f(indt,indvar) - ctdsalt_a);
  CPcor_new.sqr=CPcor_new.f(indsqr); % ds2=nanmean(s_new.f(indt,indsqr) - ctdsalt_a)
				     % Note if use ds1/2, indt is only last last from loop. 

  
  % %%%%%%% Including Cécile Cabanes' routine (from DM_CPcor) %%%%%%%%%%%%%%
  if cecile
    % -- Adapt variables to, and use, the code of COMPUTE_new_CPcor_brian: --
    CPcor_SBE=CPcor_SBE;
    CTcor_SBE=d;
    CPcor_DEF=CPcor_new.rec;
    minPRESS = min(plim); %1500; % 1000 too shallow in the Nordic Seas.
    maxTEMP = -0.5; %E
    %1. Load Argo data (e.g. 1st ascending profile if deep enough)
    %				: (1xN_LEVELS) psal_argo,  temp_argo, pres_argo, from PSAL, TEMP, PRES,with bad QCS removed (NaN values)
    %				: (1xN_LEVELS) optional : tempad_argo, presad_argo, from TEMP_ADJUSTED and PRES_ADJUSTED, if available
    %				: (1x1)        lon_argo, lat_argo  from  LONGITUDE, LATITUDE
    iia=~isnan(argo_s)&~isnan(argo_p)&~isnan(argo_t);
    psal_argo=argo_s(iia)';
    pres_argo=argo_p(iia)';
    temp_argo=argo_t(iia)';
    psal_argo_corr=argo_s(iia)';
    pres_argo_corr=argo_p(iia)';
    temp_argo_corr=argo_t(iia)';
    lon_argo=flt.lon(np);
    lat_argo=flt.lat(np);
    %2. Load reference data (eg deployment CTD)
    %				: (1xN_LEVELS2)  psal_ref, temp_ref, pres_ref,
    %				: (1x1)          lon_ref, lat_ref
    ~isnan(ctd_salt)&~isnan(ctd_z)&~isnan(ctd_temp);
    psal_ref=ctd_salt(ans)';
    temp_ref=ctd_temp(ans)';
    pres_ref=ctd_z(ans)';
    lon_ref=ctd.lon(ind);
    lat_ref=ctd.lat(ind);
    % 3. Back off temp and pressure correction values to get uncorrected float conductivity (cond_ZERO)
    a1 = (1 + CTcor_SBE.*temp_argo + CPcor_SBE.*pres_argo);   %  take raw temp & pres values
    cond_argo = gsw_C_from_SP(psal_argo,temp_argo,pres_argo);
    cond_ZERO = cond_argo.*a1;
    %4. Re-calculate salinity data by using adjusted pressure and compute Argo data’s derived quantities
    % sa : absolute salinity
    % ct : conservative temperature
    psal_argo_corr = gsw_SP_from_C(cond_argo,temp_argo_corr,pres_argo_corr);
    sa_argo_corr = gsw_SA_from_SP(psal_argo_corr,pres_argo_corr,lon_argo,lat_argo);
    ct_argo_corr = gsw_CT_from_t(sa_argo_corr,temp_argo_corr,pres_argo_corr);
    %5. Compute reference data’s derived quantities
    sa_ref = gsw_SA_from_SP(psal_ref,pres_ref,lon_ref,lat_ref);
    ct_ref = gsw_CT_from_t(sa_ref,temp_ref,pres_ref);
    %6. Compute the conductivity that the float should have used to calculate
    %   and report practical salinity that is in agreement with reference data
    %   => cond_expected
    % Interpolation of psal_ref onto float conservative temperature levels.
    [psal_ref_i,pres_ref_i] = interp_climatology(psal_ref',ct_ref',pres_ref',psal_argo_corr,ct_argo_corr,pres_argo_corr); % routine OW (=>deal with temp inv)
    % and the conductivity that the float should have used to calculate and report practical salinity that is in agreement with reference data:
    cond_expected = gsw_C_from_SP(psal_ref_i',temp_argo_corr,pres_argo_corr);
    % 7. Levels selection
    kok = find(pres_argo_corr > minPRESS & isfinite(cond_expected));
    %kok = find(temp_argo_corr < maxTEMP & isfinite(cond_expected)); %E
    % 8. Solve least square problem to get optimized values of CPcor and M (from B.King)
    % We are looking for M and CPcor so that:
    % cond_expected = (cond_ZERO * M)/(1 + t * CTcor_SBE + p * CPcor)
    % so CPcor * p - M * (cond_ZERO/cond_expected) = - (1 + t*CTcor_SBE)
    % This is a least squares problem similar to linear regression;
    % Borrow QR factorisation from polyfit;
    % In vector terms below:
    % v * [CPcor; M] = b;
    p = pres_argo_corr(kok);
    rat = -cond_ZERO(kok)./cond_expected(kok);
    b = -1 - CTcor_SBE*temp_argo_corr(kok);
    v = [p(:)/1e8 rat(:)];
    b = b(:);
    [Q R] = qr(v);
    coefs = R\(Q'*b);
    cpcor_n = coefs(1); % cpcor * 1e8
    M_n = coefs(2); % best cfac for this profile;
    CPcor_new.king = cpcor_n*1e-8;
    M_new = M_n;
    % -- From subroutine change_cpcor: --
    a1 = (1 + CTcor_SBE.*temp_argo + CPcor_SBE.*pres_argo);  %take raw temp & pres values
    cond_argo = gsw_C_from_SP(psal_argo,temp_argo,pres_argo);
    cond_ZERO = cond_argo.*a1;
    %Calculate the optimized float conductivities using  CPcor_new value and a multiplicative calibration value M_new.
    b1 = (1 + CTcor_SBE.*temp_argo_corr + CPcor_new.king.*pres_argo_corr);
    cond_argo_DM = M_new.*cond_ZERO./b1;
    % compute the corresponding psal
    psal_argo_DM = gsw_SP_from_C(cond_argo_DM,temp_argo_corr,pres_argo_corr);

    % -- Adapt to operator_CPcor_new: --
    % Put back into correct places relative to original argo data.
    [s_new.gain,pt_new.gain]=deal(nan(size(argo_s)));
    s_new.gain(iia)=psal_argo_DM;
    pt_new.gain(iia) = gsw_pt_from_t(psal_argo_DM,temp_argo_corr,pres_argo_corr,0);
    %s_new.gain=psal_argo_DM;
    %pt_new.gain = gsw_pt_from_t(s_new.gain,temp_argo_corr,pres_argo_corr,0);

  
  else
    CPcor_new.king = [];
    M_new = 1;
  end
  CPcor_new.gain=CPcor_new.king;

  
  % ------- Sal and ptemp profiles for recommended and Cécile's CPcor: --------------
  switch tool
   case 'sw'
    s_new.recl = argo_sbe_CPcor(argo_p,argo_t,argo_s,CPcor_new.rec-1.5e-8); 
    pt_new.recl = sw_ptmp(s_new.recl,argo_t,argo_p,0);
    s_new.rec = argo_sbe_CPcor(argo_p,argo_t,argo_s,CPcor_new.rec); 
    pt_new.rec = sw_ptmp(s_new.rec,argo_t,argo_p,0);
    s_new.recu = argo_sbe_CPcor(argo_p,argo_t,argo_s,CPcor_new.rec+1.5e-8); 
    pt_new.recu = sw_ptmp(s_new.recu,argo_t,argo_p,0);
    s_new.king = argo_sbe_CPcor(argo_p,argo_t,argo_s,CPcor_new.king); 
    pt_new.king = sw_ptmp(s_new.king,argo_t,argo_p,0);
   case 'gsw'
    Cnew.recl = Co.*(1 + d*argo_t + CPcor_SBE*argo_p) ./ (1 + d*argo_p + CPcor_new.rec-1.5e-8*argo_p); 
    s_new.recl = gsw_SP_from_C(Cnew.recl,argo_t,argo_p);
    pt_new.recl = gsw_pt_from_t(s_new.recl,argo_t,argo_p,0);
    Cnew.rec = Co.*(1 + d*argo_t + CPcor_SBE*argo_p) ./ (1 + d*argo_p + CPcor_new.rec*argo_p); 
    s_new.rec = gsw_SP_from_C(Cnew.rec,argo_t,argo_p);
    pt_new.rec = gsw_pt_from_t(s_new.rec,argo_t,argo_p,0);
    Cnew.recu = Co.*(1 + d*argo_t + CPcor_SBE*argo_p) ./ (1 + d*argo_p + CPcor_new.rec+1.5e-8*argo_p); 
    s_new.recu = gsw_SP_from_C(Cnew.recu,argo_t,argo_p);
    pt_new.recu = gsw_pt_from_t(s_new.recu,argo_t,argo_p,0);
  end

  

  %%%%%%%%% Figures: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  f(1)=figure(1); set(f(1),'position',figz,'name',['Compared salinity profiles ',float_names{I}]);
  %plot(argo_s,pt,'b-',ctd_salt,ctd_ptmp,'r-')
  plot(argo_s,argo_p,'g-',ctd_salt,ctd_z,'k-'); axis ij; xlabel PSAL; ylabel PRES; grid
  ylim([0 4000]);
  hleg.str=legend(leg.float,leg.ctd,'location','southwest'); title('Compared salinity profiles');
  ax(1)=gca;
  
  f(2)=figure(2); set(f(2),'position',figz,'name',['Var and sqrt ',float_names{I}]); 
  ax(5)=subplot(2,1,1); plot(CPcor_new.f*1e8,ds_var,'-',CPcor_new.f(indvar)*1e8,ds_var(indvar),'ro'); 
  xlabel('CPcor_{new} [10^{-8} dbar^{-1}]','interpreter','TeX'); ylabel var;
  ax(6)=subplot(2,1,2); plot(CPcor_new.f*1e8,ds_sqr,'-',CPcor_new.f(indsqr)*1e8,ds_sqr(indsqr),'ro'); 
  xlabel('CPcor_{new} [10^{-8} dbar^{-1}]','interpreter','TeX'); ylabel sqr;
  
  %indz=find(argo_p>=500); cindz=find(ctd_z>=500); % Only plot under 500 m
  %indz=find(pt<-0.5); cindz=find(ctd_ptmp<-0.5); % Only plot colder than -0.5C
  %indz=find(argo_p>=min(plim)); cindz=find(ctd_z>=min(plim)); % Only plot under depth corresponding to -0.5C

  f(3)=figure(3); set(f(3),'position',figz,'name',['Comparison of methods ',float_names{I}]);
  gp=plot(smoothdata(argo_s(1:end),'movmean',15,'omitnan'),pt(1:end),'g-',...
	  smoothdata(s_new.recl(1:end),'movmean',15,'omitnan'),pt_new.recl(1:end),'c--',...
	  smoothdata(s_new.rec(1:end),'movmean',15,'omitnan'),pt_new.rec(1:end),'c-',...
	  smoothdata(s_new.recu(1:end),'movmean',15,'omitnan'),pt_new.recu(1:end),'c--',...
	  smoothdata(s_new.king(1:end),'movmean',15,'omitnan'),pt_new.king(1:end),'m.',...
	  smoothdata(s_new.gain(1:end),'movmean',15,'omitnan'),pt_new.gain(1:end),'m-',...
	  smoothdata(s_new.f(1:end,indvar),'movmean',15,'omitnan'),pt_new.f(1:end,indvar),'b-',...
	  smoothdata(s_new.f(1:end,indsqr),'movmean',15,'omitnan'),pt_new.f(1:end,indsqr),'r-',...
	  smoothdata(ctd_salt(1:end),'movmean',15,'omitnan'),ctd_ptmp(1:end),'k-');%cindz
  xlabel PSAL; ylabel PTMP; title('Comparison of methods');grid
  %ylim; ylim([ans(1) -.4]); tightaxis([],'x');
  legend(gp([1 3 5:end]),...
	 ['Float PSAL / OldCPcorr (' num2str(CPcor_old*1.0e+8,'%6.2f') ' \cdot 10^{-8} dbar^{-1})'],...
	 ['NewCPcorr-recommended (' num2str(CPcor_new.rec*1.0e+8,'%6.2f') ' \cdot 10^{-8} dbar^{-1})'],...
	 ['NewCPcorr-king (' num2str(CPcor_new.king*1.0e+8,'%6.2f') ' \cdot 10^{-8} dbar^{-1})'],...
	 ['NewCPcorr-king & gain (' num2str(CPcor_new.king*1.0e+8,'%6.2f') ' \cdot 10^{-8} dbar^{-1} & ' num2str((M_new-1)*1e3,'%6.2f') ' ‰)'],...
	 ['NewCPcorr-var (' num2str(CPcor_new.var*1.0e+8,'%6.2f') ' \cdot 10^{-8} dbar^{-1})'],...
	 ['NewCPcorr-sqr (' num2str(CPcor_new.sqr*1.0e+8,'%6.2f') ' \cdot 10^{-8} dbar^{-1})'],...
	 'CTD','location','best');
  ax(2)=gca;
  
  f(4)=figure(4); set(f(4),'position',figz,'name',['Temperature and salinity ',float_names{I}]);
  ax(3)=subplot(1,2,1); gp=plot(pt(1:end),argo_p(1:end),'g-',...
				pt_new.recl(1:end),argo_p(1:end),'c--',...
				pt_new.rec(1:end),argo_p(1:end),'c-',...
				pt_new.recu(1:end),argo_p(1:end),'c--',...
				pt_new.king(1:end),argo_p(1:end),'m.',...
				pt_new.gain(1:end),argo_p(1:end),'m-',...
				pt_new.f(1:end,indvar),argo_p(1:end),'b-',...
				pt_new.f(1:end,indsqr),argo_p(1:end),'r-',...
				ctd_ptmp(1:end),ctd_z(1:end),'k-');%cindz
  axis ij; xlabel PTMP; ylabel PRES; title('Resulting potential temperatures');
  hlim=line([-.95 -.5],[1 1]*min(plim)); set(hlim,'linestyle','--','color','k');
  ylim; ylim([500 ans(2)]); tightaxis([],'x');grid
  hleg.alt=legend(gp([1 3 5:end]),...
		 ['Float PSAL / OldCPcorr (' num2str(CPcor_old*1.0e+8,'%6.2f') ' \cdot 10^{-8} dbar^{-1})'],...
		 ['NewCPcorr-recommended (' num2str(CPcor_new.rec*1.0e+8,'%6.2f') '\pm1.5 \cdot 10^{-8} dbar^{-1})'],...
		 ['NewCPcorr-king (' num2str(CPcor_new.king*1.0e+8,'%6.2f') ' \cdot 10^{-8} dbar^{-1})'],...
		 ['NewCPcorr-king & gain (' num2str(CPcor_new.king*1.0e+8,'%6.2f') ' \cdot 10^{-8} dbar^{-1} & ' num2str((M_new-1)*1e3,'%6.2f') ' ‰)'],...
		 ['NewCPcorr-var (' num2str(CPcor_new.var*1.0e+8,'%6.2f') ' \cdot 10^{-8} dbar^{-1})'],...
		 ['NewCPcorr-sqr (' num2str(CPcor_new.sqr*1.0e+8,'%6.2f') ' \cdot 10^{-8} dbar^{-1})'],...
		 'CTD','location','southeast');
  %subplot(1,2,2), plot(ctd_salt,ctd_z,'k-',s_new.f(:,indsqr),argo_p,'r-',argo_s,argo_p,'g-');
  ax(4)=subplot(1,2,2); plot(argo_s(1:end),argo_p(1:end),'g-',...
		       s_new.recl(1:end),argo_p(1:end),'c--',...
		       s_new.rec(1:end),argo_p(1:end),'c-',...
		       s_new.recu(1:end),argo_p(1:end),'c--',...
		       s_new.king(1:end),argo_p(1:end),'m.',...
		       s_new.gain(1:end),argo_p(1:end),'m-',...
		       s_new.f(1:end,indvar),argo_p(1:end),'b-',...
		       s_new.f(1:end,indsqr),argo_p(1:end),'r-',...
		       ctd_salt(1:end),ctd_z(1:end),'k-');%cindz
  axis ij; xlabel PSAL; ylabel PRES; title('Resulting salinities');
  hlim=line([34.9 34.925],[1 1]*min(plim)); set(hlim,'linestyle','--','color','k');
  ylim; ylim([500 ans(2)]); tightaxis([],'x');grid

  % The overview of figures:
  f(5)=figure(5); set(f(5),'position',scr.*[1 1 .62 1],'name',['operator CPcor_new ',float_names{I}]);
  apos=nan(4,4); 
  for i=1:4, a=subplot(2,2,i); apos(i,:)=a.Position; delete(a); end; clear a; 
  a{1}=copyobj([ax(1) hleg.str],5); set(a{1}(1),'position',apos(1,:)); title(a{1}(1),['a) ',get(a{1}(1),'title').String])
  a{2}=copyobj([ax(2)         ],5); set(a{2}(1),'position',apos(2,:)); title(a{2}(1),['b) ',get(a{2}(1),'title').String])
  a{3}=copyobj([ax(3) hleg.alt],5); set(a{3}(1),'position',apos(3,:)); title(a{3}(1),['c) ',get(a{3}(1),'title').String])
  a{4}=copyobj([ax(4)         ],5); set(a{4}(1),'position',apos(4,:)); title(a{4}(1),['d) ',get(a{4}(1),'title').String])
  a{5}=copyobj([ax(5)         ],5); set(a{5}(1),'position',[.32 .16 .13 .05]);
  a{6}=copyobj([ax(6)         ],5); set(a{6}(1),'position',[.32 .22 .13 .05]);
  set(a{5}(1),'ytick',[]); set(a{6}(1),'xtick',[],'ytick',[]); delete(get(a{6}(1),'xlabel'));
  get(a{3}(2),'position'); set(a{3}(2),'position',[.4 .466 ans(3:4)]); % The legend
  xlim([a{2}(1),a{4}(1)],[34.9 34.925]); 
  ylim([a{3}(1),a{4}(1)],[1000 3500]);
  ylim([a{2}(1)],[-.95 -.5]);
  xlim([a{3}(1)],[-.95 -.5]);
  %grid([a{1}(1),a{2}(1),a{3}(1),a{4}(1)])
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  % ------- Implement the decision (after plotting all alternatives): --
  CPcor_new=CPcor_new.(decision)

  % ------- Calculate values for for all profiles (matrix): ------------
  switch tool
   case 'sw'
    psal_new = argo_sbe_CPcor(flt.pres,flt.temp,flt.psal,CPcor_new);
    ptmp_new = sw_ptmp(psal_new,flt.temp,flt.pres,0);
    pden_new = sw_pden(psal_new,flt.temp,flt.pres,0);
   case 'gsw'
    Cnew.mat = Co.*(1 + d*flt.temp + CPcor_SBE*flt.pres) ./ (1 + d*flt.pres + CPcor_new*flt.pres); 
    psal_new = gsw_SP_from_C(Cnew.mat,flt.temp,flt.pres);
    ptmp_new = gsw_pt_from_t(psal_new,flt.temp,flt.pres,0);
    pden_new = gsw_pot_rho_t_exact(psal_new,flt.temp,flt.pres,0);
  end


  % ------- Print and save stuff: --------------------------------------
  save([cpcordir,filesep,'operatorCPcor'],'ctd_*','CPcor_new','M_new','leg','np','minPRESS');
  for i=1:4
    figure(i); print('-dpng',[cpcordir,filesep,get(gcf,'name'),'.png']);
  end
  figure(5); print(gcf,'-dpng',[outfiles{I}(1:end-4),'_CPcor_new.png']);
  
end % Float loop

  % f(5)=figure(5); set(f(5),'position',figz,'name',['Hovmoeller PSAL ',float_names{I}]);
  % [a1,a2]=argo_plott_snitt (flt.pres,psal_new,flt.juld,'psal');
  % set(a1,'CLim',[34.84 35.02])
  % f(6)=figure(6); set(f(6),'position',figz,'name',['Hovmoeller TEMP ',float_names{I}]);
  % [a1,a2]=argo_plott_snitt (flt.pres,ptmp_new,flt.juld,'temp');
  % brighten(0.5)
  % set(a1,'CLim',[-0.9 5])
  % f(7)=figure(7); set(f(7),'position',figz,'name',['Hovmoeller PDENS ',float_names{I}]);
  % [a1,a2]=argo_plott_snitt (flt.pres,pden_new-1000,flt.juld,'pden');
  % brighten(0.5)
  % set(a1,'CLim',[27.5 28.1])

  % f(8)=figure(8); set(f(8),'position',figz,'name',['Time series ',float_names{I}]);
  % for i=1:length(flt.psal(1,:))
  %   indz=find(flt.pres(:,i)>3300 & flt.pres(:,i)<3500);
  %   s_mean(i) = mean(psal_new(indz,i),'omitnan');
  %   t_mean(i) = mean(flt.temp(indz,i),'omitnan');
  % end
  % subplot(2,1,1), plot(flt.juld,s_mean); datetick
  % subplot(2,1,2), plot(flt.juld,t_mean); datetick



