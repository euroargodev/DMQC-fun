function [pres,temp,sal,long,lat,dates,dt,dmon,mons,yrs] = inpolygon_referencedata(LO,LA,tardir,maxprof,time)
% INPOLYGON_REFERENCEDATA	Finds RDB data inside lon/lat-polygon.
% 
% [pres,temp,sal,long,lat,dates,dt,dmon,mons,yrs] 
%		= inpolygon_referencedata(LO,LA,tardir,maxprof,time)
% 
% LO,LA	  = polygon within which to find all reference data.
% tardir  = the directory of referencedata, as organised by the OWC 
%	    toolbox. 
% maxprof = Integer giving the maximum number of profiles to deliver
%	    (default = all). If given, a random selection of almost
%	    this number of profiles will be given. When less profiles
%	    are found, all will be given.
%           If  maxprof < 6, then maxprof denotes the maximum monthly
%           seasonal offset allowed in the set. I.e., 0 allows for
%           reference data ±15 days around mean time, while positive
%           integers adds 30 days in both directions.
% time    = Time(s) to center the reduced selection around, in serial
%           days. When this is given, the maxprof number of profiles
%           nearest in time of year (regardless of year) will be
%           selected, in order to account for seasonality. When several
%           times are given, their mean is used. Time can also be used
%           without maxprof, to get dt output for all profiles inside
%           polygon.
%
% pres, temp, sal = matrices of chosen data.
% long, lat       = row vectors of positions of chosen data.
% dates		  = row vector of dates of chosen data (units as in
%                   reference data).
% dt		  = time differences to the central time (in days). 
% dmon            = month differences to the central time (calendar months). 
% mons		  = Reference data months of year
% yrs		  = Reference data years
%
% Note that reference data nearest in seasonality, will appear in the
% last column of the output. This is done so that the nearest will
% appear on top when plotting the output.
%
% This is part of DMQC-fun.  Requires PADCONCATENATION (provided with
% DMQC-fun), MATLAB_OWC, and EVENMAT (see INIT_DMQC).
%
% See also INPOLYGON PADCONCATENATION HALO FINDWMO

% DMQC-fun by J. Even Ø. Nilsen, Ingrid Angel, Birgit Klein, and Kjell Arne Mork.
% Main programmer: jan.even.oeie.nilsen@hi.no.

if any(contains(tardir,'historical_ctd')), tartyp={'ctd'}; end
if any(contains(tardir,'historical_bot')), tartyp(end+1)={'bottle'}; end
if any(contains(tardir,'argo_profiles')), tartyp(end+1)={'argo'}; end

error(nargchk(3,5,nargin));
if nargin<5 | isempty(time),	time=[];	end
if nargin<4 | isempty(maxprof),	maxprof=[];	end

maxprof=floor(maxprof);

wmosq=unique(findwmo(LO,LA));

[PRES,TEMP,SAL,LONG,LAT,DATES,dt,dmon,mons,yrs]=deal([]);
for j=1:length(tardir)
  for i=1:length(wmosq)
    try 
      load([tardir{j},filesep,tartyp{j},'_',int2str(wmosq(i))]);
      % for h=1:size(pres,2)
      % 	[~,IA]=unique(floor(pres(:,h)/subsample));
      % 	PRES=padconcatenation(PRES,pres(IA,h),2);
      % 	TEMP=padconcatenation(TEMP,temp(IA,h),2);
      % 	SAL=padconcatenation(SAL,sal(IA,h),2);
      % end
      PRES=padconcatenation(PRES,pres,2);
      TEMP=padconcatenation(TEMP,temp,2);
      SAL=padconcatenation(SAL,sal,2);
      LONG=[LONG,long]; 
      LAT=[LAT,lat];
      DATES=[DATES,dates];
    end
  end
end

LONG>180; LONG(ans)=LONG(ans)-360;		% Subtract 360 to match the input polygons

IN = find(inpolygon(LONG,LAT,LO,LA));		% Indices inside polygon
% Now, if there is no time input and number of found indices is less
% than maxprof, no more selection is done.

if ~isempty(IN)
  if ~isempty(time)				% PRIORITISE AROUND TIME
    rtime=rdtime(DATES(IN));		% Reference data (RD) times
    [~,~,tdays]=yrday(time(:));		% Input time in circular yeardays
    [~,~,rdays]=yrday(rtime(:));	% RD circular yeardays
    cday=nanmean(tdays);		% Central day of input
    cday=cday/abs(cday);		% Retain unit amplitude
    %[~,I]=sort(cday-rdays);		% Sort RD by distance in circle (as unique as in angle)
    % Vector based angles θ = cos-1 [ (a · b) / (|a| |b|) ]
    dt=366/2/pi*acos(dot(repmat([real(cday),imag(cday)],length(rdays),1),[real(rdays),imag(rdays)],2));%./(abs(cday).*abs(rdays)));
    [~,I]=sort(dt);			% Sort RD by angles, i.e. days differences
    IN=IN(I);				% Apply sorting to the indices for data inside polygon
    dt=dt(I);				% Apply sorting to the time differences
    rtime=rtime(I);			% Apply sorting to the RD times
    if ~isempty(maxprof)		% Reduce to the desired number:
      if maxprof==0		% Data from the same months of year
	[~,MT]=datevec(time);
	[~,MR]=datevec(rtime); 
	ismember(MR,MT);
      elseif maxprof < 6	% Maxprof denotes the maximum monthly offset
	dt<(maxprof-1)*30+15;
      else			% Just pick the maxprof nearest times
	1:min(length(IN),maxprof); 
      end
      IN=IN(ans); dt=dt(ans);
    end
  elseif ~isempty(maxprof) & length(IN)>maxprof	% JUST REDUCE TO THE DESIRED NUMBER
    IN = unique(ceil(rand(1,maxprof)*length(IN)));% by picking profiles at random
  end

  IN=flipdim(IN(:),1);		% Change order so nearest data will be plotted
  dt=flipdim(dt(:),1);		% last, i.e., on top of more offset data. 
  
  % Set the output:
  pres=PRES(:,IN);
  temp=TEMP(:,IN);
  sal=SAL(:,IN);
  long=LONG(IN);
  lat=LAT(IN);
  dates=DATES(IN);
  
  [yrs,mons]=datevec(rdtime(dates));	% Reference data years and months
  
  dt=dt(:)';				% Ensure correct orientation
  dmon=floor(dt/30.437);		% Deviation from middle in average months 
else
  [pres,temp,sal,long,lat,dates,dt,dmon,mons,yrs] = deal([]);		% Otherwise empty
end
