function [qq,nn,dd,dc]=SymStat(HAM,varargin)
% function [qq,nn,dd,dc]=SymStat(HAM [,opts])
% Wb,Sep15,16
% -------------------------------------------------------------------- %
% Examples:
% SymStat(HAM,'fout',['./SymStat_' HAM.ops(1).op.info.qtype '.dat'],'-x','-X','-V');
% -------------------------------------------------------------------- %

  getopt('INIT',varargin); o={};
     if getopt('-v'), vflag=1; o{end+1}='-v';
     elseif getopt('-V'), vflag=2; o{end+1}='-V'; else vflag=0; end

     pflag=getopt('-p');
     xflag=getopt('-x');
     Xflag=getopt('-X');
     kx   =getopt('kx',[]);

     fout=getopt('fout','');
  getopt('check_error');

  q=HAM.info.host.cgs.store;
  s=getenv('RC_STORE');
  if ~isequal(regexprep(s,'\.fs',''),q)
      wblog('NB!','changing RC_STORE');
      setRCStore(q);
  end

  L=numel(HAM.mpo); kk=2:L-1; kk(kx)=[];
  for k=kk
     HK(k)=load_dmrg_data(HAM,k,'HK');       if ~isempty(HK(k).Q)
     [QQ{k},dd{k},dc{k}]=getQDimQS(HK(k),2); end
  end

  QQ=cat(1,QQ{:});
  dd=cat(2,dd{:})';
  dc=cat(1,dc{:});

  [QQ,I,D]=uniquerows(QQ);

  [D,is]=sort(D,'descend');
  QQ=QQ(is,:); I=I(is);

  is=[I{:}];

  dd=mat2cell(dd(is),D,1);
     for i=1:numel(dd), dd{i}=mean(dd{i}); end
  dd=[dd{:}]';

  dc=mat2cell(dc(is),D,1);
     for i=1:numel(dc), dc{i}=mean(dc{i}); end
  dc=[dc{:}]';

  Hk=HK(fix(end/2));
  sym=Hk.info.qtype; nsym=numsym(Hk);

  if nargout>1, nn=D'; end

  if ~nargout, vflag=vflag+10; end

  if nsym>1 || isempty(regexp(Hk.info.qtype,'^SU\d'))
     if q, wblog('WRN',['ignoring option -p\n' ... 
          '(intended for a single SU(N) symmetry only)']);
        pflag=0;
     end
     if ~isempty(fout), wblog('WRN',['ignoring option fout\n' ... 
          '(intended for a single SU(N) symmetry only)']);
        fout='';
     end
  elseif ~nargout, pflag=1; end

  if ~isempty(fout)
    fid=fopen(fout,'w'); if fid<0, fout
       wberr('failed to open file !?'); end
    fapp=['>>' fout]; fprintf(1,'\n');

    fprintf(fid,'# %s #\n',repmat('-',1,72));
    fprintf(fid,'# L=%g, kx=%s\n',L,mat2str(kx));

    QA=uniquerows([QQ; fliplr(QQ)]);

    no=numel(HAM.ops); qloc=cell(1,no); qop=cell(1,no);
    for k=1:no
       Q=HAM.ops(k).op.Q;
          if length(Q)==2, Q{3}=zeros(size(Q{1})); end
       H0(k,:)=Q;
    end
    Sloc=cell2mat(H0); s=size(Sloc);
    Sloc=struct('Q',{mat2cell(Sloc,s(1),repmat(s(2)/3,1,3))},'qdir',[+1,-1,-1]);

    qop =uniquerows(cat(1,Sloc.Q{3}));
    qloc=uniquerows(cat(1,Sloc.Q{1:2}));

    fprintf(1,'%s\n',repmat('-',1,85));
    [AL,s1]=SymStore(sym,'--tensor-prod',QA,qloc,'--fid',fid,o{:});
    [SX,s2]=SymStore(sym,'--tensor-prod',QA,qop,'--fid',fid,o{:},'-op');

    if Xflag
    [AO,s3]=SymStore(sym,'--tensor-prod',AL.Q{2},qop,'--fid',fid,o{:},'-op');
    else s3={}; end

    ss={s1{:},s2{:},s3{:}}; ns=numel(ss);
    if ns
       fprintf(1,'\n');
       if ns<=15
            fprintf(1,'%s',ss{:});
       else fprintf(1,'%s',ss{1:5}); fprintf('  (%g entries)\n',ns-10);
            fprintf(1,'%s',ss{end-4:end});
       end
    end
    fprintf(1,'%s\n',repmat('-',1,85));

    if xflag, o{end+1}='-x'; end

    fprintf(fid,'\n#! SYM=%s\n',sym);

    print_qset(fid,AL.Q{1},'QA');
    print_qset(fid,AL.Q{2},'QX');
    print_qset(fid,AL.Q{3},'qloc');
    print_qset(fid,qop,'qop');
    fprintf(fid,'\n');

    AR=AL; p=[2 1 3];
    AR.Q=AR.Q(p); AR.qdir=AR.qdir(p);

    X4=SymStore(sym,'--contract',AL,3,Sloc,2,o{:},...
      '--fid',fid,'istr','X4=contract(AL,3,Sloc,2)');
    X3=SymStore(sym,'--contract',SX,'13*',X4,'14',o{:},...
      '--fid',fid,'istr','X3=contract(SX,''13*'',X4,''14'')');
    X4=SymStore(sym,'--contract',AR,3,Sloc,2,o{:},...
      '--fid',fid,'istr','X4=contract(AR,3,Sloc,2)');
    X3=SymStore(sym,'--contract',SX,'13*',X4,'24',o{:},...
      '--fid',fid,'istr','X3=contract(SX,''13*'',X4,''24'')');

    X4=SymStore(sym,'--contract',AL,1,SX,2,o{:},...
      '--fid',fid,'istr','X4=contract(AL,1,SX,2)');
    X3=SymStore(sym,'--contract',AL,'13*',X4,'32',o{:},...
      '--fid',fid,'istr','X3=contract(AL,''13*'',X4,''32'')');
    X4=SymStore(sym,'--contract',AR,2,SX,2,o{:},...
      '--fid',fid,'istr','X4=contract(AR,1,SX,2)');
    X3=SymStore(sym,'--contract',AR,'23*',X4,'32',o{:},...
      '--fid',fid,'istr','X3=contract(AR,''13*'',X4,''32'')');

    fprintf(1,'\n');
    return

  elseif vflag
    fprintf(1,'\n');
    ss={ '---------------------------------'
         '   symmetry #occ   dim  <dim>'
         '    labels        irep   mult'
         '---------------------------------'};

    nI=numel(I); ns=numel(ss);
    n2=ceil((nI-ns)/2);

    for k=1:n2+ns, ii=[k-ns k+n2];
       for l=1:2, i=ii(l);
          if i<=0
             ss{l,2}=ss{i+ns,1};
          elseif i<=nI
             if nsym==1
                  qs=sprintf( '%g',QQ(i,:)); qs=['(' qs ')'];
             else qs=sprintf(' %g',QQ(i,:)); qs=['[' qs ' ]'];
             end
             ss{l,2}=sprintf('%2d. %6s %4d %6g %6.1f',i,qs,D(i),dc(i),dd(i));
          else
             ss{l,2}='';
          end
       end

       fprintf(1,'   %-35s|  %-35s\n',ss{1:2,2});
    end

  end

  if pflag
    ah=smaxis(2,2,'tag',mfilename); header('%M'); addt2fig Wb
    setax(ah(1,1))

    MB=tril(ones(size(QQ,2)));

    YY=QQ*MB;

    for i=1:size(YY,2)
       [y,I,d]=uniquerows(YY(:,i));
       plot(y,d,'o-','Disp',...
          sprintf('row #%d',i));
       hold on
    end

    xlabel('length of row in Young diagram')
    ylabel('number of occurances');
    legdisp('Location','N');

  setax(ah(1,2))

    [y,I,d]=uniquerows(sum(YY,2));
    plot(y,d,'o-'); hold on

    postext('SE',['total number of boxes' 10 'in Young diagram'])
    xlabel('# boxes');
    ylabel('number of occurances');

  setax(ah(2,1))

    for i=1:size(QQ,2)
       [q,I,d]=uniquerows(QQ(:,i));
       h=plot(q,d,'o-','Disp',sprintf('q_%g',i)); hold on
       if i==1
          set(h,'MarkerS',8);
       end
    end
    set(h,'LineSt','none'); sms(h,6);

    label('value of label q_i','number of occurances');
    legdisp('Location','NE');
  end

end

% -------------------------------------------------------------------- %
function q=qmat2cellstr(q)

    if ~isnumeric(q), wberr('invalid usage'); end

    i=find(q>9); q(i)=q(i)+7; q=char(q+'0');
    q=cellstr(q)';

end

% -------------------------------------------------------------------- %

function print_qset(fid,Q,vn)

   Q=uniquerows(Q);
   i=find(Q>9); Q(i)=Q(i)+7; Q=char(Q+'0');
   qs=cellstr(Q)';

   if numel(qs)>10, s='\n#'; else s=''; end

   fprintf(fid,['\n# %s got %g symmetry sectors:' s '%s\n'],...
   vn,numel(qs),sprintf(' %s',qs{:}));
end

% -------------------------------------------------------------------- %
