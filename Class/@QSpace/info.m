function info(A,vstr,cflag)
% function info(A [,istr,cflag])
%
%    header routine used in display()
%    Options: istr (info string) and cflag (compact info).
%
% Wb,Mar01,08

  if ~nargin, helpthis, return; end

  nA=numel(A);
  nm=inputname(1); if isempty(nm) && nA>1
    if nargin>1
         nm=regexprep(vstr,'[= ]*$','');
    else nm='ans'; end
  end

  if nargin<2,
     if ~isempty(nm) && nA>1, vstr=sprintf('%s = ',nm);
     else vstr=''; end
  elseif ~ischar(vstr)
     helpthis, error('Wb:ERR','invalid usage (istr must be string)')
  end

  if nargin<3,                cflag=0;
  elseif isequal(cflag,'-C'), cflag=4;
  elseif isequal(cflag,'-c'), cflag=2;
  else cflag=1;
  end

  if cflag, nl=''; else nl=char(10); end

  for k=1:nA
     if nA>1, vstr=sprintf('%s(%g) = ',nm,k); end

     if isempty(A.Q) && isempty(A.data)
        if isempty(vstr)
             fprintf(1,[nl '  %s (empty QSpace)\n'],vstr);
        else fprintf(1,[nl '  (empty QSpace)\n']); end
     elseif cflag>1
          info_1line(A(k),vstr,cflag);
     else info_1(A(k),vstr,cflag); end
  end

end

% -------------------------------------------------------------------- %
function info_1line(A,vstr,cflag)

  r=length(A.Q); nd=numel(A.data); sx={};

  if isempty(vstr), s0='';
  else s0=sprintf('%-6s ',vstr); end

  if ~isfield(A.info,'qtype') || isempty(A.info.qtype)
       sym='(U1)';
  else sym=['''' A.info.qtype '''']; end

  stags=''; ocflag=0;
  ftags=sprintf('%%-%gs',4+4*r);

  if isfield(A.info,'itags') && ~isempty(A.info.itags)
     q=getqdir(A); ocflag=all(q>0);
     stags=itags2str(A.info.itags);
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
     sdc=sprintf([Dfmt ' %s'],r,sdc);
     sdim='';
  else
     sdc=class(A.data{1});
     if zflag && isequal(sdc,'double'), sdc='complex'; zflag=0; end
     sdc=sprintf([Dfmt ' %s'],r,sdc);

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

  s=sprintf(['%s %-6s %-12s ' ftags ' %8s  %10s  %s'],...
    s0,sym,sdc,stags,snrm,sbytes,sdim);

  if ocflag && exist('printfc')==3
       printfc(['\e[38;5;12m' s '\e[0m\n']);
  else fprintf(1,'%s\n',s);
  end

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
     s{end+1}=[ itags2str(A.info.itags) ];
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

