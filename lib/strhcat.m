function ss=strhcat(varargin)
% Function: s=strhcat(ss [,sep])
%
%    Concatenate strings horizontally with
%    given separator (default: blank)
%
% Options
%
%   '-s',sep    alternative way of specifying separator (blank)
%   '-a'        show strings as argument list 'var1', 'var2', 'var2', ...
%
% Wb,Aug01,07  Wb,Apr08,08

% NB! see also matlab's strjoin(c, ', ') introduced in matlab/2013a
% Starting in matlab/2016, the `join' function with is also intended
% for strings anyway with similar syntax! join(c,', ')

  if ~nargin, ss=''; return; end
  if nargin==2 && iscell(varargin{1})
       sep=varargin{2}; varargin=varargin(1);
  else sep={}; end

  getopt('init',varargin);
     sep  =getopt({'sep','-s'},sep);
     aflag=getopt('-a');  % print as argument list: 'var1', 'var2', 'var2',...
  varargin=getopt('get_remaining'); narg=length(varargin);

  if aflag
     if ~isempty(sep)
          sep=['''' sep ''''];
     else sep=''', '''; end
  elseif isempty(sep), sep=' '; end

  if narg
     if iscell(varargin{1}), ss=varargin{1}; e=(narg>1);
     else e=0; ss=varargin;
        for i=1:length(ss)
        if ~ischar(ss{i}), e=1; break; end, end
     end
  else e=1; end

  if e
     eval(['help ' mfilename]);
     if nargin || nargout, wberr('invalid usage'), end, return
  end

  if iscell(ss)
       if isempty(ss), ss=''; return; end
       s=strtrim(cat(2,ss{:}));
  else s=strtrim(reshape(ss,1,[]));
  end
  if isempty(s), ss=''; return; end

  if ischar(ss)
     ss(:,end+1)=10;
     ss=strread(ss','%s','whitespace',' \n');
  end

  ss=reshape(ss,1,[]);
  ss(2,:)={sep}; ss{end}='';

  ss=cat(2,ss{:}); if aflag, ss=['''' ss '''']; end

end

