function [dd,qq,I]=getSymDim(varargin)
% function [dd,qq,I]=getSymDim(varargin)
%    
%    usage: getSymDim <copy/paste string of CData/CGRef info>
% => use regexp to extract symmetry type and labels of input
%
% Wb,Jan16,15

  for i=1:nargin, if ~ischar(varargin{i})
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end, end

  t=varargin; t(2,1:end-1)={' '}; istr=[t{:}];
  t=strread(istr,'%s','whitespace',' ')';

  sym=regexprep(t{1},'(S[pU])\(*(\d+)\)*','$1$2');

  t=t(2:end); if ~isempty(t)
     t(2,1:end-1)={' '}; istr=[t{:}];
     istr=regexprep(istr,'_\d+','');
     qs=strread(istr,'%s','whitespace',' (),;*');
     for i=1:numel(qs)
        if isempty(regexp(qs{i},'^[\d\w]+$'))
        wbdie('invalid input string'); end
     end
     nq=numel(qs);
  else nq=0; end

  D=get_dir(sym,'RStore');

  rsym=get_rank(sym); vflag=0;

  xb=(1:rsym)';

  if ~nq
     ff=dir(D);
       ff=ff(find([ff.isdir]==0));
       [~,is]=sort([ff.bytes]); ff=ff(is);
     nq=numel(ff); qs=cell(1,nq); vflag=2;

     for i=1:nq
        qs{i}=regexprep(ff(i).name,'\(([\w\d\s]+)\).*','$1');
     end
  else
     ff=cell(1,nq);
     for i=1:nq, ff{i}=dir([D '/*' qs{i} '*']); end
     ff=[ff{:}]; nq=numel(ff);
  end

  if ~nargout && vflag<=0, vflag=1; end
  if vflag, fprintf(1,'\n'); end

  istr1=''; istr2=''; istr3=''; istr={''}; dmax1=1;
  err=zeros(1,nq);
  er1=zeros(1,32);

  for i=1:nq
     mat=sprintf('%s/(%s).rep',D,qs{i});
     if ~exist(mat,'file'), wbdie('file ''%s'' not found',mat); end
     I=load(mat,'-mat'); q=I.RSet;
     if isfield(I,'err')
        err(i)=I.err;
        l=find(q.J); if numel(l)==1, er1(l,q.J(l))=I.err; end
     end

     d=size(q.Z,1); t=regexprep(q.type,'[\(\)]','');
     if vflag
        if all(q.J>=0 & q.J==round(q.J))
             qstr=qmat2char(q.J);
        else qstr=vec2str(q.J,'-f'); end

        istr={sprintf('%s  (%s)  d=%d',t,qstr,d)};

        if vflag>1, istr{2}=datestr(ff(i).datenum);
               fprintf(1,'  %3d.  %-24s   %s\n',i,istr{:});
        elseif fprintf(1,'  %3d.  %s\n',i,istr{1});
        end
     else qstr='';
     end

     nb=q.J*xb;
     if i==1 || nb>nbmax, nbmax=nb; istr2=istr{1}; end
     if i==1 || d>dmax, dmax=d; istr3=istr{1}; end

     S.qstr=qstr;
     S.qset=q.J;
     S.dim=d;
     S.nY=nb;
     S.istr=istr{1};

     IR(i)=S;
     if i==1, IR(nq).qset=[]; end
  end

  dd=cat(1,IR.dim);
  qq=cat(1,IR.qset);

  if isempty(t)
     [dd,is]=sort(dd); qq=qq(is,:); IR=IR(is);
  end

  if vflag && ~isempty(regexpi(sym,'SU')), fprintf(1,'\n');

     nb=[IR.nY]'; N=rsym+1;
     br=mod(nb,N);
     bs=sparse(br+1,ones(size(br)),ones(size(br)),N,1);

     fprintf(1,'# statistics in number of boxes in young tableau mod N=%g ',N);
     fprintf(1,'(max nY=%g)\n',max(nb));
     fprintf(1,'# n_boxes%%%g   n_irep | histogram\n',N); s1={};
     for i=1:N
        s1{i}=sprintf(' %6g  %6g',i-1,full(bs(i)));
     end
     s1=cat(1,s1{:});

     n=max(nb)+1; e=ones(length(nb),1);
     bs=full(sparse(nb+1,e,e,n,1))'; n=length(bs); m=min(8,max(bs)); q=m/max(bs);
     s2=repmat(' ',n,m);
     for i=1:n
        l=ceil(bs(i)*q); s2(i,1:l)='*';
     end
     s=repmat(' ',n,1); i=1:10:n; s(i)='0'+(i-1)/10;
     s2=fliplr([s,s2])'; s2=[repmat(' ',size(s2,1),2) s2];

     l=[size(s1,1),size(s2,1)];
     if l(1)<l(2)
        n=size(s1,2);
        s1=[ repmat(' ',1,n); s1; repmat(' ',diff(l)-1,n) ];
     elseif l(1)>l(2)
        n=size(s2,2);
        s2=[ repmat(' ',diff(l),n), s2 ];
     else fprintf(1,'\n');
     end

     disp([s1,repmat(' ',size(s1,1),6),s2]);
     fprintf(1,'\n');

     if 1
        M=ones(size(qq,2)); M=M-triu(M,1);
        nn=qq*M;

        n1=nn(1,:); rsym=length(n1); k=0; dx=zeros(rsym,0);

        if norm(n1)==0
           while 1, d=zeros(size(n1)); k=k+1; n1_=n1;
              for j=1:rsym
                 if j==1 || n1(j)<n1(j-1), n1(j)=n1(j)+1;
                    l=matchIndex(nn,n1);
                    if ~isempty(l)
                       d(j)=dd(l);
                       dx(j,n1(j))=d(j);
                    else n1(j)=n1(j)-1; end
                 end
              end

              dmax=max([d,dmax]);
              if all(d==0), break; end
           end
           l=matchIndex(nn,n1);
           q1=n1*inv(M);
           istr1=IR(l).istr;
        else
           n1=[]; q1=[]; istr1='';
        end
     end

     if ~isempty(istr1)
        s='got all Young tableauxs within';
        fprintf(1,'# %-46s %s (d<=%g)\n',s,istr1,dmax1);
     end

     if ~isempty(istr2)
        s=sprintf('Young tableaux with largest number of boxes',nbmax);
        fprintf(1,'# %-46s %s\n',s,istr2);
     end
     if ~isempty(istr3)
        s='largest irrep in RStore';
        fprintf(1,'# %-46s %s\n',s,istr3);
     end

     i=find(err); e=err(i);
     if ~isempty(e)
        s='range of CR accuracy';
        fprintf(1,'\n# %-46s %.3g .. %.3g (%d/%d)\n',s,min(e),max(e),numel(e),nq);
     end

     q=sum(er1~=0,1); l=find(q<size(er1,1),1); j=find(q);
     i=find(er1); e=er1(i);
     if ~isempty(e)
        s=sprintf('... for box tableauxs (%d..%d cols)',l-1,max(j));
        fprintf(1,'# %-46s %.3g .. %.3g (%d/%d)\n',s,min(e),max(e),numel(e),nq);
     end

     qq=cat(1,IR.qset); r=size(qq,2);
     fprintf(1,'\n# %-46s (i=1..%g)\n','label range for each q_i',r);
     for j=1:size(qq,2)
        fprintf(1,'       q%d  %s\n',j,vec2str(unique(qq(:,j)),'-c'));
     end

     fprintf(1,'\n');
  elseif ~nargout, fprintf(1,'\n');
  end

  if ~nargout, clear dd
  elseif nargout>2
     if vflag
          I=add2struct('-',IR,q1,n1,dx);
     else I=IR; end
  end

end

