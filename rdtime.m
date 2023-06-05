function t = rdtime(dates)
% RDTIME	Translates time format in Argo reference data
%		to Matlab serial days.
%
% t = rdtime(dates)
% 
% dates = the time variable in Argo reference database.
%
% t     = same times in Matlab serial days.
%
% The time in serial days is a direct integer version of a datetime
% string (ISO 8601) 'yyyymmddTHHMMSS', i.e. DATESTR format 30, without
% 'T'.
% 
% See also DATENUM2 DATESTR

% Last updated: Wed May 24 11:33:38 2023 by jan.even.oeie.nilsen@hi.no

% 30 (ISO 8601)   'yyyymmddTHHMMSS'       20000301T154517 
% Original datenum only takes formats 0,1,2,6,13,14,15,16,23

D=size(dates);

if ~isempty(dates)
  d=int2str(dates(:));				% Convert integer to string dates
  t=datenum2(strcat(d(:,1:8),'T',d(:,9:14)),30);% Add the T and convert
  t=reshape(t,D);
else
  t=[];
end

