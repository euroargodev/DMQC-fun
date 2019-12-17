function wmo = findwmo(lon,lat)
% FINDWMO	Gives the WMO square numbers of positions.
% 
% wmo = findwmo(lon,lat)
% 
% Returns NaNs for NaN positions. 
% 
% See also CREATE_WMO_BOXES

error(nargchk(2,2,nargin));
size(lon);
if ans~=size(lat), error('Input is not of equal size!'); end

wmo=nans(ans);

lon>180;lon(ans)=lon(ans)-360;
if any(lon<-180) | any(180<lon) | any(lat<-90) | any(90<lat)
  error('Non-existent position entered!');
end

wmo(lon>=0 & lat>=0) = 1000;  % 1=NE
wmo(lon>=0 & lat< 0) = 3000;  % 3=SE
wmo(lon< 0 & lat< 0) = 5000;  % 5=SW 
wmo(lon< 0 & lat>=0) = 7000;  % 7=NW

lon(lon==180)=179; lat(lat==90)=89;

floor(lat/10); ans(ans<0)=ans(ans<0)+1; wmo=wmo+abs(ans)*100;
floor(lon/10); ans(ans<0)=ans(ans<0)+1; wmo=wmo+abs(ans);

