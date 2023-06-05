% PAIR_FLOATS_WITH_REFERENCEDATA compares float profiles to individual
% reference data profiles in the vicinity of float.
%
% DMQC-fun v0.9.
% J. Even Ø. Nilsen, Ingrid M. Angel-Benavides, Birgit Klein, Malgorzata Merchel, and Kjell Arne Mork.
% Last updated: Wed May 24 15:27:57 2023 by jan.even.oeie.nilsen@hi.no

% PAIR_FLOATS_WITH_REFERENCEDATA is mainly for DMQC in shallow seas. It
% finds referencedata near in time and space to any float positions,
% visualises and calculates simple salinity offsets that can be used in
% judging acccuracy of the float salinity in areas where there are no
% deep waters and OWC calibration is not possible.
%
% The calculated offsets are only shown in the figures, they are not
% stored for any other use. The only action that can come out of these
% judgements is to decide that part of the time series is bad. This is
% done by setting cal_action=[0 4] for this float in INIT_DMQC and
% chosing which cycles is bad by editing calseries in the local
% set_calseries.m.
%
% For the report this script makes two sections, one section to
% replace the OWC subsection in the 'Correction of Salinity Data'
% section, and one appendix with all the plots. If OWC has been
% completed, the former will not be made (or erased by
% RUN_OW_CALIBRATION). 
%
% This function operates on all floats selected in INIT_DMQC.
%
% Ideally wait until all floats are prepared, so that data for
% comparisons become valid, but will also work with unfinished visual
% checking.
%
% Set your mdx and mdt in the init section below!
% 
% ----------------------------------------------------------------------

clear all; close all
init_dmqc; % Paths and filenames etc. for chosen set of floats.

% ----- Init: ----------------------------------------------------------
mdx=15;			% max distance allowed, in km
mdt=15;			% max time difference allowed, in days
% ----------------------------------------------------------------------

o=length(float_names);	% Number of floats

for K=1:o		% Loop floats

  % Go to that float's working dir:
  float_working_dir=[my_working_dir,'DMQC',filesep,float_names{K}];
  cd(float_working_dir);	
  
  % Load float data:
  load(outfiles{K},'LAT','LONG','DATES','PRES','SAL','TEMP','PTMP','PROFILE_NO','JULD','CYCLE_NUMBER');
  LONG>180; LONG(ans)=LONG(ans)-360;
  T=datenum(DATES,1,1);

  % Find all reference data in the area the float traverses:
  find(~isnan(LONG)&~isnan(LAT));
  [LO,LA]=halo(LONG(ans),LAT(ans),.6,.4);
  [pres,temp,sal,long,lat,dates]=inpolygon_referencedata(LO,LA,tardir);
  t=rdtime(dates);

  [m,n]=size(PRES);
  
  % 3D-figure of positions:
  nearpos(LONG,LAT,T,long,lat,t,mdx,mdt);
  set(gcf,'name',['pointcloud_',float_names{K}]);
  
  % Find groups of nearest neighbours:
  [dx,dt,I,J,JJ] = nearpos(LONG,LAT,T,long,lat,t,mdx,mdt);
  I=unique(I);		% Float profiles with coincidences
  M=length(I);		% Number of float profiles with coincidences

  
  % ----- The section in the report: -----------------------------------
  
  if M==0 % No pairings found
    fid=fopen('coincidence.tex','w'); % Empty the file, i.e., no section in report
    fida=fopen('coincidence_appendix.tex','w'); % Empty the file, i.e., no appendix in report
  else
    % Make new section file for the report (if no OWC done):
    fido=fopen('owc.tex','r'); s=fscanf(fido,'%s'); fclose(fido);
    if startsWith(s,'Figure')
      fid=fopen('coincidence.tex','w');
      % Title and text:
      fprintf(fid,'%s\n','\subsubsection{Comparison with near-coinciding reference data}\label{sec:coincidence}');
      fprintf(fid,'%s%s%s%s%s%s%s%s%s%s\n',['Float ',float_names{K},' has no OWC made, instead single float profiles have been compared to reference data profiles that happens to be in the vicinity. For Float ',float_names{K},' there are ',int2str(M),' cycles with available reference data within ',int2str(mdx),'~km and ', int2str(mdt),'~days: ',zipnumstr(I),' (see Appendix~\ref{app:coincidence}']);  
      if M==1, fprintf(fid,'%s\n',[' Figure~\ref{fig:coincidence',int2str(M),'}).']);
      else,    fprintf(fid,'%s\n',[' Figures~\ref{fig:coincidence',int2str(1),'}--\ref{fig:coincidence',int2str(M),'}).']);
      end
      fprintf(fid,'%s\n',' ');
      fprintf(fid,'%s\n',['The purpose of the comparisons are not so much the calculated', ...
			  ' salinity offsets, as that would require individual judgement', ...
			  ' as to whether there are comparable layers. Shallow waters are', ...
			  ' quite often not layered stably and quite often has strong', ...
			  ' variability in space and time, so reference profiles not exactly', ...
			  ' coinciding with a float profile will not be relevant. On the', ...
			  ' other hand in cases were both profiles have stable layers and', ...
			  ' shared potential-temperature levels in those', ...
			  ' layers, calculation or judgement of offset in those layers', ...
			  ' may be relevant.']);
    end % section for report    
    
    % Make new appendix file for the report (also when OWC is done, if any pair):
    fida=fopen('coincidence_appendix.tex','w');
    % Title and text:
    fprintf(fida,'%s\n','\subsection{Comparison with near-coinciding reference data}\label{app:coincidence}');
    fprintf(fida,'%s\n',['For Float ',float_names{K},' there are ',int2str(M),' cycles with available reference data within ',int2str(mdx),'~km and ', int2str(mdt),'~days: ',zipnumstr(I),' ']);  
    if M==1, fprintf(fida,'%s\n',[' (Figure~\ref{fig:coincidence',int2str(M),'}).']);
    else,    fprintf(fida,'%s\n',[' (Figures~\ref{fig:coincidence',int2str(1),'}--\ref{fig:coincidence',int2str(M),'}).']);
    end
    fprintf(fida,'%s\n',' ');
    fprintf(fida,'%s\n',['These figures are merely meant to provide some guidance to', ...
			' whether a float may have drifting salinities or other obvious errors.', ...
			' No quantification can be gained from these results, without carefully', ...
			' performed comparisons on a case by case basis.']);

    % Make the figure files and floats in report:
    for i=1:M		% Loop float profiles with coincidence

      % Start figure for float cycle:
      fprintf(fida,'%s\n','\begin{figure}[H]');
      fprintf(fida,'%s\n','\centering');
      
      for j=1:min(length(JJ{I(i)}),3) % Loop 3 nearest reference profiles 
      
	% Display message:
	fltncyc=['Float ',float_names{K},' cycle ',int2str(CYCLE_NUMBER(I(i)))];
	disp(['Finding simple offset for ',fltncyc,': ', ...
	      num2ordinal(j),' nearest reference profile found ',int2str(round(dx(I(i),JJ{I(i)}(j)))),...
	      ' km away and ',int2str(round(dt(I(i),JJ{I(i)}(j)))),' days past.']);
	
	% Calculate salinity offset (also makes the figure files):
	[sal_offset] = salinity_offset(PRES(:,I(i)),TEMP(:,I(i)),SAL(:,I(i)),pres(:,JJ{I(i)}(j)),temp(:,JJ{I(i)}(j)),sal(:,JJ{I(i)}(j)),[5 inf],[-2 5],LONG(I(i)),LAT(I(i)),T(I(i)),long(JJ{I(i)}(j)),lat(JJ{I(i)}(j)),t(JJ{I(i)}(j)),fltncyc,'CTD');
	
	% Print produced figure to file:
	filetag=['salinity_offset_','cycle',int2str(CYCLE_NUMBER(I(i))),'neighbour',int2str(j)];
	%print('-dpng',[outfiles{K}(1:end-4),'_',filetag,'.png'])
	print('-depsc',[outfiles{K}(1:end-4),'_',filetag,'.eps'])
      
	% Add figure to section for report:
	fprintf(fida,'%s%s%s\n','\includegraphics[width=\textwidth,natwidth=1200,natheight=500]{\floatsource\WMOnum_',filetag,'.eps}');
	
      end % neighbour loop
	
      % End figure for float cycle:
      fprintf(fida,'%s\n',['\caption{Float \WMOnum\ Cycle~',int2str(CYCLE_NUMBER(I(i))),' and near coinciding reference data profile(s). ',...
			  'Left panel(s): Salinity profiles against pressure. ',...
			  'Right panel(s): Salinity profiles against potential temperature.}']);
      fprintf(fida,'%s%s%s\n','\label{fig:coincidence',int2str(i),'}');
      fprintf(fida,'%s\n','\end{figure}');
    
    end   % Loop float profiles with coincidence
  end   % If any float profiles with coincidence
  fclose(fid); 
  fclose(fida); 
end % loop floats

cd(my_working_dir); % Go back to the main working dir
