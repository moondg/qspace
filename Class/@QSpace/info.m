function info(A,varargin)
% function info(A [,istr,cflag,rmax])
%
%    short general info on QSpace
%    e.g., used as QSpace header in display().
%
% Options
%
%    istr   info string
%    cflag  whether to show info in compact form ('-c','-C')
%    rmax   max. rank to consider (used to reserve space for itags
%           so that entries line up over repeated callls e.g.
%           for a QSpace array.
%
% Wb,Mar01,08

  getopt('INIT',varargin);
     if     getopt('-c'), cflag=2;
     elseif getopt('-C'), cflag=4;
     else   cflag=getopt('cflag',0); end

     if getopt('~oc'), ocflag=0; else ocflag=1; end

  args=getopt('get_remaining'); nargs=numel(args);

  if nargs && isnumeric(args{end})
       rmax=args{end}; nargs=nargs-1;
  else rmax=0; end

  vstr='';
  if nargs==1, vstr=args{1};
     if ~ischar(vstr) || ~isempty(vstr) && vstr(1)=='-'
        disp(vstr); wbdie('invalid istr');
     end
  elseif nargs
     if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  nm=inputname(1); nA=numel(A);
  if isempty(nm) && nA>1
     if ~isempty(vstr)
          nm=regexprep(vstr,'[= ]*$','');
     else nm='ans'; end
  elseif isempty(vstr) && nA>1
     vstr=sprintf('%s = ',nm);
  end

  if cflag, nl=''; else nl=char(10); end

  for k=1:nA
     if nA>1, vstr=sprintf('%s(%g) = ',nm,k); end

     if isempty(A.Q) && isempty(A.data)
        if isempty(vstr)
             fprintf(1,[nl '  %s (empty QSpace)\n'],vstr);
        else fprintf(1,[nl '  (empty QSpace)\n']); end
     elseif cflag>1
          info_1line(A(k),vstr,cflag,ocflag,rmax);
     else info_1(A(k),vstr,cflag); end
  end

end

% -------------------------------------------------------------------- %
function info_1line(A,vstr,cflag,ocflag,rmax)

  rk=length(A.Q); nd=numel(A.data); sx={};

  if isempty(vstr), s0='';
  else s0=sprintf('%-6s ',vstr); end

  if ~isfield(A.info,'qtype') || isempty(A.info.qtype)
       sym='(U1)';
  else sym=['''' A.info.qtype '''']; end

  stags='';

  if isfield(A.info,'itags') && ~isempty(A.info.itags)
     if ocflag, q=getqdir(A); ocflag=all(q>0); end
     stags=itags_to_str(A.info.itags,'QS:info');
  end

  if isfield(A.info,'otype') && ~isempty(A.info.otype)
     sx{end+1}=A.info.otype;
  end

  s=A.data; s=whos('s'); sbytes=num2str2(s.bytes,'-b');

  if cflag, xsep='x'; Dfmt='%dD'; else xsep=' x '; Dfmt='%d-D'; end

  zflag=~isreal(A);

  if ~nd
     sdc=class(A.data);
     if zflag && isequal(sdc,'double'), sdc='complex'; zflag=0; end
     sdc=sprintf([Dfmt ' %s'],rk,sdc);
     sdim='';
  else
     sdc=class(A.data{1});
     if zflag && isequal(sdc,'double'), sdc='complex'; zflag=0; end
     sdc=sprintf([Dfmt ' %s'],rk,sdc);

     if ~isempty(A.Q) && ~isempty(A.Q{1})
        dd=getDimQS(A);
        if isvector(dd)
           dstr=vec2str(dd,'sep',xsep,'fmt','%d','nofac','-f');
        else
           dstr=vec2str(dd(1,:),'sep',xsep,'fmt','%d','nofac','-f');
           for i=2:size(dd,1)
               if ~isequal(dd(i,:),dd(i-1,:))
               q=int2str2(dd(i,:)); dstr=[dstr, ' | ', strhcat(q,xsep)]; end
           end
        end
     else dstr=''; end

     sdim=dstr;
  end

  if zflag
     if cflag, sx{end+1}='C';
     else sx{end+1}='complex'; end
  end

  if numel(A.data)==1 && numel(A.data{1})==1
       a=A.data{1};
  else a=normQS(A); end
  snrm=num2str(a,'%.4g');

  if ~isempty(sx)
     if cflag
          sdim=sprintf('%-16s%s',sdim,sprintf(' %s',sx{:}));
     else sdim=sprintf('%-16s%s',sdim,sprintf(' %10s',sx{:}));
     end
  end

  e1=''; em=''; 
  if wblog('--hl-check')
     if ocflag
        e1=[char(27) '[38;5;12m'];
        em=[char(27) '[0m'];
     end

     l=4+5*max(3,rmax);
     q=regexprep(stags,'\x1B\[[\d;]+m','');
     q=diff([length(q), l]);
     if q>0, stags = [stags, repmat(' ',1,q)]; end
  else
     l=4+6*max(3,rmax);
     q=diff([length(stags), l]);
     if q>0, stags = [stags, repmat(' ',1,q)]; end
  end

  fprintf(1,[ e1 '%s %-6s ' em stags  e1 ' %-12s %8s  %10s  %s' em '\n'],...
  s0,sym,sdc,snrm,sbytes,sdim);

end

% -------------------------------------------------------------------- %

function info_1(A,s,cflag)

  if ~isempty(s), fprintf(1,'\n%s\n',s);
  elseif ~cflag,  fprintf(1,'\n'); end

  scalar=0;

  if isempty(A.Q)
     if numel(A.data)==1 && numel(A.data{1})==1, scalar=1; end
     fprintf(1,'     Q:  []');
  else
     s=size(A.Q{1});
     for i=2:numel(A.Q)
        if ~isequal(s,size(A.Q{i})), wbdie('invalid Q data'); end
     end
     fprintf(1,'     Q:  %dx {%2d x %d }',numel(A.Q),s);
  end

  s={};
  if isfield(A.info,'otype') && ~isempty(A.info.otype)
     s{end+1}=A.info.otype;
  end
  if isfield(A.info,'itags') && ~isempty(A.info.itags)
     s{end+1}=[ itags_to_str(A.info.itags,'QS:info') ];
  end

  if ~isreal(A), s{end+1}='complex'; 
      if scalar, s{end}=[s{end} ' scalar']; end
  elseif scalar, s{end+1}=[ class(A.data{1}), ' scalar'];
  end

  if ~isempty(s), s(2,:)={',  '}; s=s([2 1],:); s{1}(1)=[]; end
  s=cat(2,s{:});

  if ~isfield(A.info,'qtype') || isempty(A.info.qtype), s={'',s};
     if ~isempty(A.Q) && size(A.Q{1},2)>1
          s{1}='all abelian U(1)';
     else s{1}='abelian U(1)'; end
     fprintf(1,'  %s%s\n',s{:});
  else
     fprintf(1,'  having ''%s''%s\n',A.info.qtype,s);
  end

  r=length(A.Q);
  s=A.data; s=whos('s'); s=num2str2(s.bytes,'-b');

  s={s,''};

  if ~cflag && isfield(A.info,'ctime')
     if getenvb('QS_DEBUG')
        [cts,q]=ctime(A.info.ctime,'-s');
        if abs(q)>1, s{2}=['    ' cts]; end
     end
  end

  if cflag
       xsep='x';   Dfmt='%dD';
  else xsep=' x '; Dfmt='%d-D'; end

  if scalar
     fprintf(1,['  data:  { %s }\n\n'], num2str(A.data{1}));
  elseif isempty(A.data) || isempty(A.data{1})
     fprintf(1,['  data:  ' Dfmt ' %s (%s) %s\n\n'], r, class(A.data), s{:});
  else
     if ~isempty(A.Q) && ~isempty(A.Q{1})
        dd=getDimQS(A);
        if isvector(dd)
           dstr=vec2str(dd,'sep',xsep,'fmt','%d','nofac','-f');
        else
           dstr=vec2str(dd(1,:),'sep',xsep,'fmt','%d','nofac','-f');
           for i=2:size(dd,1)
               q=int2str2(dd(i,:)); dstr=[dstr, ' => ', strhcat(q,xsep)]; 
           end
           q=normQS(A); q(2)=q*q;
           if q<0.01, i=[]; else i=find(abs(q-round(q))<1E-8); end
           if ~isempty(i), q=rat2(q(1)); else q=sprintf('%.4g',q(1)); end
           dstr=[dstr sprintf('  @ norm = %s',q)];
        end
     else dstr=''; end

     fprintf(1,['  data:  ' Dfmt ' %s (%s)      %-16s%s\n'], ...
       r, class(A.data{1}), s{1}, dstr, s{2});
     if ~cflag, fprintf(1,'\n'); end
  end

end

% -------------------------------------------------------------------- %

