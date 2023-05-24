function h=postext(varargin)
% function h=postext(xn,yn,fmt,.. [,{topts}]))
% function h=postext({pos},fmt,.. [,{topts}]))
%
%    Add sprintf formatted text in the current axis
%    at position (xn,yn) in normalized units
%    with (0,0) the lower left corner.
%    where fmt, .. is forwarded to sprintf()
%
%    Instead of xn and yn, the position may also be specified
%    similar as in legdisp() using postrans().
%
%    The last argument, if it is a cell, is assumed
%    to contain text options applied to the text handle h.
%
% Wb 2000, Wb,Mar31,07

  tag='postext';
  if nargin==1 && regexp(varargin{1},'^-+rm')
     delete(findall(gca,'tag',tag));
     return
  end

  if nargin<2, eval(['help ' mfilename]); return; end

  o={};

  if isnumeric(varargin{1})
     if isnumeric(varargin{2})
          x=varargin{1};    y=varargin{2};    varargin(1:2)=[];
     else x=varargin{1}(1); y=varargin{1}(2); varargin(1  )=[];
     end
  else
     [pos,o]=postrans(varargin{1}); varargin(1)=[];
     x=pos(1); y=pos(2);
  end

  getopt('init',varargin);
     rmflag=getopt('-rm');
     dx    =getopt('dx',0);
     totop =getopt('-top');
  args=getopt('get_remaining'); nargs=length(args);

  if isempty(args), eval(['help ' mfilename]); return; end
  if ~isempty(dx)
     if length(dx)>2, wberr('invalid dx'); end
     x=x+dx(1); if length(dx)>1, y=y+dx(2);
  end

  if rmflag
     delete(findall(gca,'type','text','tag',tag));
     if ~nargs, return; end
  end

  topts={};
  if nargs>1 && iscell(args{end})
     topts=args{end}; args(end)=[];
     nargs=numel(args);
  end

  if nargs>1
     fmt=regexprep(args{1},'\\+','\\\\');
     fmt=regexprep(fmt,'\\\\n','\\n');
     s=sprintf(fmt,args{2:end});
  elseif iscell(args{1})
     s=strhcat(args{1},char(10));
  else
     s=args{1};
  end
  s=regexprep(s,'\\+','\\');

  s=regexprep(s,'\{\s*\\+n\s*\}',char(10));
  s=regexprep(s,'\\+n\>',char(10));
  s=regexprep(s,'\\+t\>',char(9));

  h=text(x,y,s,'Units','normalized','tag',tag,o{:});
  if ~isempty(topts), set(h,topts{:}); end

  if totop
     set(h,'Units','Data'); p=get(h,'Pos');
     p(3)=0.1; set(h,'Pos',p);
  end

  if ~nargout, clear h; end

end

