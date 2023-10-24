function t = rdtime(dates,opt)
% RDTIME	Translates time format in Argo reference data base
%		to Matlab serial days.
%
% t = rdtime(dates,opt)
% 
% dates = the time variable in Argo reference database (RDB).
% opt	= option string containing:
%		'reverse' : translate serial days to RDB-dates
%
% t     = same times in Matlab serial days (or RDB-dates when reversing).
%
% The time in serial days is a direct integer version of a datetime
% string (ISO 8601) 'yyyymmddTHHMMSS', i.e. DATESTR format 30, without
% 'T'.
% 
% See also DATENUM2 DATESTR

% Last updated: Fri Oct 20 17:10:05 2023 by jan.even.oeie.nilsen@hi.no

% 30 (ISO 8601)   'yyyymmddTHHMMSS'       20000301T154517 
% Original datenum only takes formats 0,1,2,6,13,14,15,16,23

error(nargchk(1,2,nargin));
if nargin < 2 | isempty(opt), opt=''; end

if ~isempty(dates)
  if contains(opt,'reverse')
    t=datestr(dates,30);	% Make ISO8601 from serial day 
    t=replace(string(t),'T','');% Remove the Ts
    t=str2num(char(t)); % Convert to float
  else
    D=size(dates);
    d=int2str(dates(:));			% Convert integer to string dates
    t=datenum2(strcat(d(:,1:8),'T',d(:,9:14)),30);% Add the T and convert
    t=reshape(t,D);
  end
else
  t=[];
end
