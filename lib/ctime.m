function [tt,dtnow]=ctime(tt,sflag)
% function [tt,dtnow]=ctime(tt [,'-s'])
%
%    convert C-time or linux system time `date +%s` to matlab format.
%    dtnow return time difference to now in seconds.
%
% Options
%
%   '-s'  convert result to string format
%
% Wb,Oct06,14

  persistent tref

  if ~nargin
     helpthis, if nargin || nargout
     wbdie('invalid usage'), end, return
  end

  if isempty(tref)
     tref=datetime('01-Jan-1970','TimeZone','local');
     tref=datenum(tref) + days(tzoffset(tref));
  end

  if nargin<2, sflag=0;
  elseif isequal(sflag,'-s'), sflag=1; 
     tn=mod(tt,1);
  else wbdie('invalid usage'); end

  tt = tt/(24*3600) + tref;

  if nargout>1
     q=now-tt; a=abs(q); i=find(a==max(a),1);
     dtnow=q(i)*24*3600;
  end

  if sflag, tt=matcell(tt);
     for i=1:numel(tt)
       if tn(i)
          tt{i}=datestr(tt{i},'mm/dd/yy HH:MM:SS.FFF');
       else
          tt{i}=datestr(tt{i},'mm/dd/yy HH:MM:SS');
       end
     end
     if numel(tt)==1, tt=tt{1}; end
  end

end

