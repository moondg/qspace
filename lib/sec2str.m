function ss = sec2str(tt,varargin)
% function ss = sec2str(tsec [,opts])
% Wb,Jul10,03

  if nargin<1 || ~isnumeric(tt)
     helpthis, if nargin || nargout
     wbdie('invalid usage'), end, return
  end

  getopt('init',varargin);
     mflag=getopt('-m');
     hflag=getopt('-h');
     dflag=getopt('-d');
  getopt('check_error');

  nt=numel(tt); ss=cell(1,nt);

  for i=1:numel(tt), t=tt(i);

     if t<60
        if t==round(t)
             ss{i}=sprintf('%g secs',t);
        else ss{i}=sprintf('%.1f secs',t);
        end
     else
        t=round(t); s=mod(t,60); m=(t-s)/60;

        if m<60 || mflag
           ss{i}=sprintf('%d:%02d',m,s);
        else
           t=m; m=mod(t,60); h=(t-m)/60;

           if dflag
              if h<10
                   ss{i}=sprintf('%.3f hrs',tt(i)/3600);
              else ss{i}=sprintf('%.2f days',tt(i)/(24*3600));
              end
           elseif h<24 || hflag
              ss{i}=sprintf('%d:%02d:%02d',h,m,s);
           else
              t=h; h=mod(t,24); d=(t-h)/24;
              ss{i}=sprintf('%g-%02d:%02d:%02d',d,h,m,s);
           end
        end
      end
   end

   if nt==1, ss=ss{1}; end

end

