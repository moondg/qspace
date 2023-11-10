function cto(varargin)
% Function cto - change to directory
% Usage: cto <pid>
% 
%    switches directory where pid is path identifier (acronym)
%    for a full directory. This looks for a perl script cto.pl.
%    If that does not exist, this looks for an environmental variable
%    with <pid> in upper case letters that defines a path.
%    If none succeed, this returns an error.
%    
% Options
%    
%    -q   quiet
%    -p   toggle to PROJECT directory
%    
% See also unix shell script: cto <dkey>
% Wb,Jun 2005  Wb,May17,06

  if ~nargin
   % eval(['help ' mfilename]);
     cd ~ ; return
  end

  getopt('init',varargin);
     vflag=~getopt({'-q','nolog'});
     pflag= getopt('-p');
  varargin=getopt('get_remaining'); narg=numel(varargin);
  
  if narg~=1 || ~ischar(varargin{1})
     helpthis; wbdie('invalid usage'); 
  elseif isequal(varargin{1},'.'), return; end

  pid=varargin{1}; p=''; p0=pwd;

  try
     [e,s]=system('command -v cto.pl');
     if ~e
      % user defined perl script that translates pid to full path
        p=evalc(['! cto.pl -d ' pid '  2>/dev/null']);
        cd(p);
     else
      % look for paths defined via environmental variables
        p=getenv(upper(varargin{1}));
        if ~isempty(p), cd(p);
        else wbdie('invalid tag ''%s''',varargin{1}); end
     end

  catch me
     if ~isempty(find(p==10))
     fprintf(1,  '   WRN directory contains NEWLINE !?'); end
     fprintf(1,'\n   ERR failed to change directory to `%s''\n\n',varargin{1});
     rethrow(me);
  end
  
  if vflag && ~isequal(p0,pwd)
   % show new directory if it actually changed
     printf('>> %s\n', pwd);
  end

end

