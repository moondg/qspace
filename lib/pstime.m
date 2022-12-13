function [s,e]=pstime(varargin)
% function [s,e]=pstime([opts])
% Wb,Mar30,19

  getopt('init',varargin);
     vflag=getopt('-v');
     sflag=getopt('-s');
  getopt('check_error');

  p=getpid();

  try
    if sflag
       [s,e]=evalc(sprintf('system(''pidstat.pl -ml %g'')',p));
       s=eval(s);
    elseif vflag
       [s,e]=evalc(sprintf('system(''pstime.pl -Vml %g'')',p));
    else
       [s,e]=evalc(sprintf('system(''pidstat.pl %g'')',p));
    end
  catch l
     wblog('ERR',''); s
     l=lasterror; dispstack(l.stack);
  end

end

