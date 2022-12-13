function setRCStore(varargin)
% function setRCStore(varargin)
% Wb,Aug29,16

  getopt('init',varargin);
     vflag=getopt('-v');
     if getopt('--leo'),
       [i,RCS]=system('cto.pl -d rcx');
     elseif getopt('--loc'),
       [i,RCS]=system('cto.pl -d rcs');
     else
       helpthis, if nargin || nargout
       error('Wb:ERR','invalid usage'), end, return
     end
  getopt('check_error');

  if i || ~exist(RCS,'dir')
     error('Wb:ERR','\n   ERR invalid RCS directory'); end

  q=getenv('RC_STORE');
  if ~isequal(q,RCS)
     fprintf(1,...
        '\n   Changing RC_STORE ...\n       %s\n   --> %s\n\n',q,RCS);
     setenv('RC_STORE',RCS);
  elseif vflag
     fprintf(1,'\n   Keeping RC_STORE = %s\n\n',RCS);
  end

end

