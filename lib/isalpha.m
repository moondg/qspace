function s=isalpha(varargin)
% Function: s=isalpha([OPTS,] s)
% Options:
%
%    '-l'  lenient (also accepts +-.= as alpha)
%
% Wb,Aug06,07

  getopt('init',varargin);
     strict=~getopt('-l');
  varargin=getopt('get_remaining');

  if length(varargin)~=1 || ~ischar(varargin{1})
     eval(['help ' mfilename]);
     if nargin || nargout, wbdie('invalid usage'), end, return
  end

  if strict
       pat=    '[0-9a-zA-Z]';
  else pat='[+=.-0-9a-zA-Z]'; end

  s=regexprep(varargin{1},pat,char(1));
  s=(double(s)==1);

end

