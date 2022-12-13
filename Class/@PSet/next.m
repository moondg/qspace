function i=next(P,varargin)
% function i=next(P [,'-v|V'])
%
%    reset(P);
%    while next(P,'-V'), structexp(P.info.p);
%       if P.info.mark(end)
%        % last index changed => got full cycle in P.s(1:end-1)
%          save(...)
%       end
%       ...
%    end
%
% Wb,Aug04,08 ; Wb,Mar13,19

  if nargin>1
     getopt('INIT',varargin);
        if getopt('-v'), vflag=1;
        elseif getopt('-V'), vflag=2; end
     getopt('check_error');
  else vflag=0;
  end

  vn=inputname(1);
  if isempty(vn), error('Wb:ERR',...
    '\n   ERR invalid usage (PSet not specified simply by name)'); end

  if P.i>=P.n, i=0;
     if P.i>P.n, wblog('WRN',...
        'reached end of PSet ''%s'' (%d/%d)',inputname(1),P.i,P.n); end
     P.i=P.i+1; assignin('caller',vn,P);
     return
  end

  if P.i<1, P.i=0; P.t=[]; end
  P.t(end+1,:)=[P.i, now];

  i=P.i+1; P.i=i;

  I=[]; ixl=[];
  if ~isempty(P.info)
     if numel(P.info)==1 && isstruct(P.info), I=P.info; ixl=I.ix;
     else wblog('WRN','overwriting %s.info data',vn); end
  end

  if isempty(I)
     I=struct('ix',[],'pstr',[],'tstart','','tfin',[],'tstr',[],'p',[],'mark',[]);
  end

  tt=P.t;
  m=prod(P.s); if m>16
     s=P.s; s=[s(1), prod(s(2:end))];
     if s(1)>2 && s(2)>8, m=s(2);
     else m=ceil(m/16); end
     if m<8, m=min(prod(P.s),8); end
  end
  tt=tt(max(1,end-m+1):end,:); nt=size(tt,1);

  if nt>1 && all(diff(tt(:,1))==1)
     np=1;

     if nt>2
        [pp,s,mu]=polyfit(tt(:,1),tt(:,2),np);
        [t,dt]=polyval(pp,P.n,s,mu);
     else
        t=tt(1,2) + (diff(tt(:,2))/diff(tt(:,1))) * (P.n-tt(1,1));
        dt=Inf;
     end

     if dt<1, dt=dt*24;
        if dt>1.4, s=sprintf(' (@ %.1f hrs',dt);
        elseif dt>0.0233
             s=sprintf(' (@ %.1f min',dt*60);
        else s=sprintf(' (@ %.1f sec',dt*3600); end
        s=[s sprintf('; %s to go)',sec2str(max(0,(t-now)*24*3600)))];
     else s=''; end

     I.tfin=t;
     I.tstr=sprintf('time finished ~ %s%s',...
       datestr(t,'ddd mm/dd/yy HH:MM:SS'),s);
  end

  [I.p,I.pstr,I.ix]=subsref(P,struct('type','()','subs',{[]}));

  if isempty(ixl)
       I.mark=ones(size(I.ix));
  else I.mark=(ixl~=I.ix);
  end

  P.info=I;

  assignin('caller',vn,P);

  if vflag && P.n>1
     if vflag>1
          banner(  '%s\n%s',I.pstr,I.tstr);
     else banner(2,'%s\n%s',I.pstr,I.tstr); end
  end

end

