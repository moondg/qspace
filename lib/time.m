function tstr = time(varargin)
% Function: tstr = time([tstamp, flag])
% default tstamp is now().
%
%    return the current date and time in the format
%    Sep10,07 14:18:38
%
% Flags
%
%    '-a'   plain datestr(now): 10-Sep-2007 14:17:12
%    '-r'   Sep 10, 2007
%    '-t'   14:17:12
%
% Wb,Nov17,99, Wb,Sep10,07

% a   = clock;        % 6dim vector ([year month day hour minute seconds])
% str = date;         % 09-Oct-2000
% str = datestr(now); % 09-Oct-2000 17:47:14

% 2002:10:08 13:42:53
% tstr = sprintf ('%d:%02d:%02d %02d:%02d:%02d', ...
%         a(1), a(2), a(3), a(4), a(5), round(a(6)) );

  if nargin && isnumeric(varargin{1})
  t=varargin{1}; varargin(1)=[]; else t=now; end
  narg=length(varargin);

  if narg==0, iflag=''; else iflag=varargin{1}; end

  switch lower(iflag)
    case '-a',  tstr = datestr(t);
    case '-r',  tstr = datestr(t,'mmmdd,yy');
    case '-l',  tstr = datestr(t,'mmm dd, yyyy HH:MM:SS');
    case {'t','-t'} ,  tstr = datestr(t,'HH:MM:SS'); % `t' for compatibility
    case '-dt', tstr = datestr(t,'yymmdd_HHMMSS');
    case '-wb', tstr = ['Wb' datestr(t,'yymmdd')];
    case '',    tstr = datestr(t,'mmmdd,yy HH:MM:SS');
    otherwise,  wbdie('invalid tag %s',iflag);
  end

end

