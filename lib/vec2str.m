function sout=vec2str(v, varargin)
% function sout=vec2str(v [,fmt,opts])
%
%    Write numerical vector as string.
%    If input is matrix, this calls mat2str2() instead.
%
% Options
%
%    'fmt',...  number format (%g)
%    'sep',..   separation of values (' ')
%    'xsep',..  customized extra set of separators ({})
%               expected in cell format: { i1,'sep1'; i2,'sep2'; ...}
%    '-f'       show full vector without any shortcuts such as x1:dx:x2
%
% Compact format notation
%
%    fmt may specify left/right bracket, as well as separator
%    in addition to format, assuming the following convention
%    '(left-boundary)(actual-fmt)(separator)(right-boundary)'
%    where space in the left-boundary is assumed replicated
%    with the right boundary.
%
%  Examples:
%      vec2str(a,'a=[ %g,  ]','-f')
%
% see also vsprintf.m, strjoin(..,delim), sizestr(), sizestrd()
% Wb,Nov12,99, Jul08,03, Jan16,06; Wb,Dec07,20; Wb,Dec07,22

  getopt ('init', varargin);

     fmt   = getopt('fmt','');
     sep   = getopt('sep','');
     rowsep= getopt('rowsep','; ');
     xsep  = getopt('xsep',{});
     vflag = getopt('-v');
     fflag = getopt('-f'); % formerly 'full' // Wb,Jul07,11
     usefac=~getopt('nofac');

  if isempty(fmt)
       fmt=getopt('get_last',fmt);
  else getopt('check_error'); end

  sout=''; % sout='[]'; // inconsistent with no brackets below!

  if issparse(v), v=full(v);
  elseif isinteger(v)
     v=double(v);
  end
  n=numel(v); 

  bL=''; bR='';
  if isempty(fmt)
     if isempty(sep), sep=' '; end
     fmt='%.5g';
  else
     fmt=regexprep(fmt,'^([^%]+)(?@bL=$1;)','');
     fmt=regexprep(fmt,'([^\w]+)$(?@bR=$1;)',''); nR=length(bR);

     if isempty(sep), sep=' ';
     elseif nR
        wblog('WRN','sep=''%s'' specified may be overwritten by fmt',sep);
     end

     if ~isempty(bL), q=''; i=0;
        regexp(bL,'[\(\[\{\\''\"=](\s+)(?@q=$1;)%');
        if ~isempty(q)
           regexp(bR,['(' q '[^\s].*)(?@i=length($1);)']);
           if i, i=nR-i; if i>1
              sep=bR(1:i); bR=bR(i+1:end); end
           end
        end
        if i<=1
           if ~isempty(regexp(bR,'[^\s]\s+([^\s].*)(?@i=length($1);)'))
              i=nR-i;
              sep=bR(1:i); bR=bR(i+1:end);
           elseif ~isempty(regexp(bR,'^(\s*[,;]\s*)(?@i=length($1);)'))
              sep=bR(1:i); bR=bR(i+1:end);
           end
        end
     else
        i=regexp(bR,'[\}\)\]]');
        if ~isempty(i),      sep=bR(1:i-1); bR=bR(i:end);
        elseif ~isempty(bR), sep=bR; bR='';
        end
     end
  end

  if n<=2
     if n, q=cell(1,n);
        for i=1:n, q{i}=sprintf_1(v(i),fmt,vflag); end
        q(2,1:end-1)={sep}; sout=[bL q{:} bR];
     end
     return
  end

  if ~isvector(v)
     o={'fmt',fmt,'sep',sep}; if fflag, o{end+1}='-f'; end
     if ~isempty(rowsep), o={o{:},'rowsep',rowsep}; end
     sout=[bL mat2str2(v,o{:}) bR];
     return
  end

  if ~fflag && usefac
     q=norm(diff(v)); if q, q=q/norm(v); end
     if q<1E-8
        sout=[ bL, sprintf_1(v(1),fmt,vflag) ...
               sprintf(' (x%g)',numel(v)), bR ];
        return
     end
  end

  if ~fflag
     if n==3 && norm(v-round(v)), fflag=2;
     elseif ~isreal(v), fflag=3; end
  end

  dxconst_flag=0;

  if ~fflag && (n>4 || ~vflag)
     dV=diff(v); ddv=std(dV);
     dv=mean(dV); adv=abs(dv);

     if ddv==0
        sout={ sprintf_1(v(1),fmt,vflag), sprintf_1(v(end),fmt,vflag) };
        if abs(dv-1)<1E-12
             sout=[bL sout{1} ':' sout{2} bR];
        else sout=[bL sout{1} ':' sprintf_1(dV(1),fmt,vflag) ':' sout{2} bR];
        end
        return
     end

     eps=1E-4*max(abs(v));
     if ddv/(adv+eps) < eps
        dxconst_flag = 1;
     end
  end

  fac=mean(abs(v)); if fac, fac=floor(log10(fac)); end
  if usefac && abs(fac)>3
       fac=10^fac; v=v/fac;
  else fac=1; end

  if dxconst_flag

     delta=(max(v)-min(v)) / (n-1);
     if delta~=0
        q=10^(round(log10(delta))+5);
        delta=(delta+q)-q;
     end

     sout={ sprintf_1(min(v),fmt,vflag), sprintf_1(max(v),fmt,vflag) };
     if abs(dv-1)<1E-12
          sout=[bL sout{1} ':' sout{2} bR];
     else sout=[bL sout{1} ':' sprintf_1(delta,fmt,vflag) ':' sout{2} bR];
     end
  else
     n=numel(v); for i=n:-1:1, sout{i}=sprintf_1(v(i),fmt,vflag); end
     sout(2,1:end-1)={sep};
     for i=1:size(xsep,1)
        sout(2,xsep{i,1})={xsep{i,2}};
     end

     sout=[bL sout{:} bR];
  end

  if fac~=1
     if ~isempty(bL)
          sout=sprintf('%G * %s%s%s',fac,bL,sout,bR);
     else sout=sprintf('%G * [%s]',fac,sout);
     end
  end

end

% -------------------------------------------------------------------- %
function s=sprintf_1(z,fmt,vflag)

  if numel(z)~=1, wberr('invalid usage'); end
  if isreal(z), s=sprintf(fmt,z); return; end

  r=real(z); i=imag(z);
  if i && abs(r/i)<1E-8, s=sprintf([fmt 'i'],i); return; end
  if r && abs(i/r)<1E-8, s=sprintf( fmt,     r); return; end

  s={ sprintf(fmt,r), sprintf(fmt,i) };

  l=[ length(s{1}), length(s{2}) ];
  w=[ s{1}([1,end])==' '; s{2}([1,end])==' ' ];
  if any(w(:))
     for i=1:2
        if s{i}(1  )==' ', s{i}=regexprep(s{i},'^\s+',''); end
        if s{i}(end)==' ', s{i}=regexprep(s{i},'\s+$',''); end
     end
  end

  if vflag
     s{2}=regexprep(s{2},'^([\+\-]?)',' $1 ');
     if s{2}(1)~=' ', s{2}=[ ' + ' s{2}]; end
  elseif isempty(regexp(s{2},'^([\+\-]?)'))
     s{2}=[ '+' s{2}];
  end

  s=[s{:} 'i']; n=min(l); l=length(s);

  if l<n
     if     any(w(:,1)), s=[repmat(' ',1,n-l) s];
     elseif any(w(:,2)), s=[s, repmat(' ',1,n-l)]; end
  end

end

% -------------------------------------------------------------------- %

