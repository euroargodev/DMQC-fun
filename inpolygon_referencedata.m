function [pres,temp,sal,long,lat,dates] = inpolygon_referencedata(LO,LA,tardir)
% INPOLYGON_REFERENCEDATA	Finds RDB data inside lon/lat-polygon.
% 
% [pres,temp,sal,long,lat,dates] = inpolygon_referencedata(LO,LA,tardir)
% 
% LO,LA	 = polygon within which to find all reference data.
% tardir = the directory of referencedata, as organised by the OWC 
%	   toolbox. 
% 
% This is part of DMQC-fun.
% Requires PADCONCATENATION (provided with DMQC-fun),
% MATLAB_OWC, and EVENMAT (see INIT_DMQC).
%
% See also INPOLYGON PADCONCATENATION HALO FINDWMO

% DMQC-fun by J. Even Ø. Nilsen, Ingrid Angel, Birgit Klein, and Kjell Arne Mork.
% v0.9.3, jan.even.oeie.nilsen@hi.no.

tartyp={'ctd_','bottle_','argo_'};

wmosq=unique(findwmo(LO,LA));

[PRES,TEMP,SAL,LONG,LAT,DATES]=deal([]);
for j=1:length(tardir)
  for i=1:length(wmosq)
    try 
      load([tardir{j},filesep,tartyp{j},int2str(wmosq(i))]);
      PRES=padconcatenation(PRES,pres,2);
      TEMP=padconcatenation(TEMP,temp,2);
      SAL=padconcatenation(SAL,sal,2);
      LONG=[LONG,long]; 
      LAT=[LAT,lat];
      DATES=[DATES,dates];
    end
  end
end

LONG>180; LONG(ans)=LONG(ans)-360;   % Subtract 360 to match the input polygons

IN = inpolygon(LONG,LAT,LO,LA);

pres=PRES(:,IN);
temp=TEMP(:,IN);
sal=SAL(:,IN);
long=LONG(IN);
lat=LAT(IN);
dates=DATES(IN);
