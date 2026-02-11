function [q,ss,qq]=LoadCData(sym,varargin)
% function [q,ss]=LoadCData(sym [,opts])
%    
%    load specified CData structure.
%
%    If the second argument specifies qdir,
%    then all CData that match given pattern are loaded.
%
% Options
%
%    '-f','-a' enforces to load full CData (and not just RefInit())
%    '-F'      same, but do *not* convert MPFR -> double
%              full size-specification (i.e., including trailing singletons)
%              can still be found in q.cgd CRef
%
%    '-x'      check whether file exists, then returns [i,full_path_to_cdata]
%    '--exist' check whether file exists, then returns [i,full_path_to_cdata]
%
% Wb,Jan30,15

  cflag=0;

  if isfield(sym,'type'), C=sym; sym=C.type; cflag=1;
  elseif nargin<2
     if nargin, sym=check_sym(sym); end ; helpthis
     if nargin || nargout, wbdie('invalid usage'); end
     return
  end

  [D0,sym]=get_dir(sym,'CStore');

  if iscell(D0), D0=D0{1}; end

  getopt('INIT',varargin);
     if getopt('-x'); xflag=1;
     elseif getopt('--exist'); xflag=2;
     else xflag=0; end

     if     getopt('-a'), aflag=2;
     elseif getopt('-f'), aflag=3;
     elseif getopt('-F'), aflag=1;
     else aflag=0; end

     if     getopt('-v'), vflag=1;
     elseif getopt('-V'), vflag=2; else vflag=0; end

     if ~nargout, vflag=vflag+10; end

  varargin=getopt('get_remaining');

  if cflag
     if ~isempty(varargin), varargin, wbdie('invalid usage'); end
     if vflag || nargout<=1, o={'-v'}; else o={}; end
     [f,Iq]=QSet2file(C,o{:}); f=[D0 '/' f];
     if ~exist(f,'file')
        wbdie('invalid CData (file not found)\n%s',f); end
     if nargout>1, ss=Iq; end

     q=load_CData_1(f,aflag);

     if vflag, s='';
        if bitand(vflag,2),
           if isequal(C.cid,q.cid)
              fprintf(1,'\n   input exactly same as loaded\n');
           elseif isequal(C.cid(1:3),q.cid(1:3))
              fprintf(1,'\n   input same as loaded (except for ID: %s)\n\n',...
              cgd_type(C.cid(4)));
           else
              fprintf(1,'\n');
              print_qset('   input',C,'-v');
           end
           s=sprintf('   loaded ( %s )',repHome(f));
        end
        print_qset(s,q,'-v');
     end
     return
  end

  if ~isempty(varargin) && ~isempty(regexp(varargin{1},'^\+*\-*$'))
     pflag=1; qdir=varargin{1}; varargin(1)=[];

     D=[D0 '/' qdir];
     if ~exist(D,'dir')
          wbdie('got invalid/non-existent directory %s',D);
     else cd(D); end

  else pflag=0;
  end

  qstr=[varargin{:}]; ff=[];

  if pflag
     if isempty(qstr), pflag=2;
     else
        ff=dir(qstr);
        ff(find([ff.isdir]))=[];
     end
  end

  if isempty(ff)
     qstr=regexprep(qstr,'.*(','');
     qstr=regexprep(qstr,').*','');
     qstr=regexprep(qstr,',(\d+)\*',';$1','once');
     qstr=regexprep(qstr,'\*','');

     if pflag, ff=dir(['*' qstr '*.cgd']); end
  end

  vlog=0;

  if pflag, nf=numel(ff); clear q
     if nf==0
        if vflag, wblog(' * ',...
          'no CData found that matches %s/*%s*',qdir,qstr); end
        return
     elseif vflag || nf>100
        s=sprintf('found %g CData that match %s/*%s*',nf,qdir,qstr);
        if nf>100, wblog('WRN','%s',s); else wblog(' * ','%s',s); end
     end
     for i=nf:-1:1
        if vflag || nf>10 || aflag
           fprintf(1,'\r... %5g/%g loading %s/%s%30s\r',i,nf,qdir,ff(i).name,'');
           vlog=vlog+1;
        end
        try   qi=load_CData_1(ff(i).name,aflag);
        catch l
           disp(l.message); if vlog, fprintf(1,'\n'); end
        end
        qi.file=ff(i);
        q(i)=qi;
     end
  else
     qdir=get_qdir(qstr);

     f=[D0 '/' qdir '/(' qstr ').cgd'];
     q=exist(f,'file'); if xflag>1, ss=f; return, end
     if q
        q=load_CData_1(f,aflag);
        q.file=dir(f);
        if xflag, ss=f; return; end

     else
        q2=split(qstr,';'); estr='';
        for i=1:numel(q2)
           qi=textscan(q2{i},'%s','whitespace',' ,()'); qi=qi{1};
           if numel(qi)>1
              [q2{i},is]=sort(qi);
              if any(diff(is)~=1)
                 estr='(input qlabels not sorted)'; break
              end
           end
        end
        wbdie('%s\nfile not found %s',repHome(f),estr);
     end
  end

  if vlog, fprintf(1,'\r%80s\r',''); end

  if nargout>1 || pflag>1, qd=[q.cgd];
     ss=double(cat2(1,qd.S,{1}));
     qq=double(cat (1,q.qset)); r=length(q(1).qset)/length(q(1).qdir);
  end

  if pflag>1 && isequal(q(1).qdir,'++--'), s2=ss;
     OM=1; if size(ss,2)>4, OM=max(s2(:,end)); else s2(:,5)=1; end

     for p=1:2, j=2*(p-1)*r+1:r:2*p*r+1;
        q2=[qq(:,j(1):j(2)-1); qq(:,j(2):j(3)-1)];
        [q1,I,D]=uniquerows(q2); [~,is]=sort(D','descend'); m=10;
        if p==1, s={'in','++'}; else s={'out','--'}; end
        fprintf(1,'\n   Most common %s multiplets (%s):\n',s{:});
        for l=1:m, i=is(l);
           fprintf(1,'%8s %d\n',QSet2str(q1(i,:)),D(i));
        end
        fprintf(1,'\n');
     end

     s2(:,1:2)=sort(s2(:,1:2),2);
     s2(:,3:4)=sort(s2(:,3:4),2);

     ah=smaxis(2,2,'tag',mfilename);
     header('%M :: sym=%s/CStore/%s',sym,qdir); addt2fig Wb
     setax(ah(1,1))

        e1=ones(size(s2,1),1);
        om=full(sparse(s2(:,end),e1,e1,OM,1));
        semilogy(om,'o'); sms(4); xtight(1.1);
        label('OM','# occurances',sprintf('max(OM)=%g',OM))

     for p=1:2
     setax(ah(2,p)), j1=2*(p-1)+1; j2=2*(p-1)+2;
        for m=1:OM
            i=find(s2(:,end)==m);
            loglog(s2(i,j1),s2(i,j2),'o','MarkerS',3+m/2); hold on
        end
        label(sprintf('s_%g',j1),sprintf('s_%g',j2))
        postext('SE','[marker size indicates OM (=3+om/2)]',{'FontSize',10});
     end

  end

  if vflag
     if pflag && nf>1
        print_qset(q,'-c');
     else
        print_qset(q,'-v');
        if isfield(q,'file')
           fprintf(1,'   file      %s\n',repHome(f));
           fprintf(1,'   size      %s\n\n',fsize2str(q.file.bytes));
        end
     end
  end

end

% -------------------------------------------------------------------- %
% outsourced // Wb,Jul04,22

function q=load_CData_1(f,aflag)
  if aflag
     I=load2(f,'-mat'); q=I.CRef; q.cdata=I.cdata;
     if aflag>1
        q.CData=cgs2double(q.cdata);
     end
  else
     I=load2(f,'CRef','-mat'); q=I.CRef;
  end
end

% -------------------------------------------------------------------- %

