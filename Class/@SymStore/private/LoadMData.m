function [M,qex]=LoadMData(sym,varargin)
% function [M,qex]=LoadMData(sym [,'-a'], varargin)
%    
%    load specified map3 data structure.
%
% Wb,Sep09,15

  D0=get_dir(sym,'-c','CStore'); qdir='++-'; qex=1;

  getopt('init',varargin);
     aflag=getopt('-a');

     if     getopt('-v'), vflag=2; 
     elseif getopt('-q'), vflag=0; else vflag=1; end

  varargin=getopt('get_remaining');

  if numel(varargin)==0
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  if ischar(varargin{1})
     q=regexprep([varargin{:}],'.*(',''); q=regexprep(q,').*','');
     q=regexprep(q,',(\d+)\*',';$1','once');
     q=regexprep(q,'\*','');
  else
     if numel(varargin)~=2 || ...
       ~isnumeric(varargin{1}) || ~isnumeric(varargin{2})
        error('Wb:ERR','\n   ERR invalid usage');
     end
     q=cat(1,varargin{:});
     q=qmat2char(q);
     q=sprintf('%s,%s',q(1,:),q(2,:));
  end

  D=append2cstr(D0,['/' qdir]); nD=numel(D);
  for l=1:nD, if ~exist(D{l},'dir')
     fprintf(1,'  WRN %s() skipping non-existing ''%s''\n',...
     mfilename,repHome(D{l}));
  end

  for l=1:numel(D)
     f=[ D{l} '/(' q ').mp3'];
     if exist(f,'file')
        M=load2(f,'-mat');
        if ~aflag
           e=sum(M.c2eps(:));
           if e>1E-30 || vflag>1 || ~nargout
              fprintf(1,'\n   Got mp3 data at e=%.3g\n\n',e);
           end
           M=M.map3;
           M.c2eps=e;
        end
        return
     end
  end

  if aflag, r=get_rank(sym);
     ff=get_files(D,['*' q '*.mp3']); nf=numel(ff);
     if nf && numel(q)==r
        w2=zeros(1,nf); d2=zeros(1,nf); j2=cell(1,nf);
        for i=nf:-1:1
           M(i)=load(ff(i).name,'-mat');

           j12=QSet2str(M(i).map3.J12);

           if isequal(q,j12(1:r))
                w2(i)=2; d2(i)=M(i).map3.cgr(1).size(2); j2{i}=j12(r+1:end);
           else w2(i)=1; d2(i)=M(i).map3.cgr(1).size(1); j2{i}=j12(1:r);
           end
        end

        [~,is]=sort(d2); i2=is(find(w2(is)==2)); n2=numel(i2);

        fprintf(1,'\n   Found %g (%g) entries ...\n',n2,nf);

        for i=1:n2, k=i2(i);
           if w2(k)==2
              fprintf(' %8g.  (%s)%8d\n',i,j2{k},d2(k));
           end
        end

        j2=uniquerows(cat(1,j2{:}));

        [jj,dd]=LoadRData(sym);

        [i1,i2,I]=matchIndex(double(jj),double(j2));
        if ~isempty(I.ix2)
           error('Wb:ERR','\n   ERR failed to find R-data in mp3 !?'); 
        end

        [~,is]=sort(dd(I.ix1)); nx=numel(is);

        fprintf(1,'\n   Entries not yet in tensorprod with ''%s'' ...\n',q);

        n1=numel(find(dd(I.ix1)<max(d2)));
        n2=max(0,min(nx-n1,5));
        for i=1:n1, k=I.ix1(is(i));
           fprintf(' %8g.  (%s)%8d\n',i,jj(k,:),dd(k));
        end
        if nx>n1+n2, fprintf(1,'\n      (skipping %g entries)\n\n',nx-n1-n2); end
        for i=nx-n2+1:nx, k=I.ix1(is(i));
           fprintf(' %8g.  (%s)%8d\n',i,jj(k,:),dd(k));
        end
     end
  elseif ~vflag, qex=0;
     q2=sort(strread(q,'%s','whitespace',',')); q=[q2{1},',',q2{2}];
     q2=char2qmat(cat(1,q2{:}));
     M=struct('J12',reshape(q2',1,[]),'J',[],'omult',[],'cgr',[]);

     fc=get_files(D,['*' q '*']);
     nf=numel(fc); q3=zeros(nf,size(q2,2)); om=zeros(nf,1);
     ss=cell(nf,1);

     for i=nf:-1:1
        I=load(fc(i).name,'CRef','-mat'); I=I.CRef; cgr(i)=I;
        q3(i,:)=I.qset(end*2/3+1:end);
        ss{i}=I.cgd.S; ss{i}(end+1:4)=1;
     end

     if nf
        ss=cat(1,ss{:});
        M.cgr=cgr; M.J=q3; M.omult=ss(:,4);
     end
  else
     error('Wb:ERR','\n   ERR file not found "%s"',f);
  end

end

