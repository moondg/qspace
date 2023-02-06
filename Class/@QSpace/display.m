function display(A,varargin)
% Function display(A [,OPTS])
% Options
%
%   '-a','-f'  show all entries (like full, force)
%   '-v'       verbose flag (shows CGS dimension separately for each symmetry)
%   '-E'       sort QSpace wrt. data (energy) [assumes diagonal Hamiltonian]
%   '-s'       sort QIDX
%   'sperm',.. sort QIDX using given permutation
%
%    By default, QSpaces with many entries (length > 14) will show
%    a truncated list, including also the two largest entries.
%
% Wb,Mar01,08

  getopt('init',varargin);
     m  = getopt('m',[6 2 2]);
     nm = getopt('nm','');

     if getopt('-v'), vflag=2;
     elseif getopt('-c'), vflag=0; else vflag=1; end

     aflag=getopt('-a'); if ~aflag,
     aflag=getopt('-f'); end
     if aflag, m=Inf; end

     if     getopt('-E'), Eflag=1; os={'-E'};
     elseif getopt('-R'), Eflag=2; os={'-R'};
     elseif getopt('-S'), Eflag=3; os={'-S'};
     else Eflag=0;
         sflag=getopt('-s');
         sperm=getopt('sperm',{});
     end

  n=getopt('get_last','');

  if isempty(nm), nm=n;
  elseif ~isempty(n), {nm,n}
     error('Wb:ERR','invalid usage (name specified twice !?)')
  end

  if Eflag, [A,isd]=sort(A,os{:}); Eflag=isd;
  else
     if ~isempty(sperm)
        if  length(A(1).Q)~=length(sperm), sperm
            error('Wb:ERR','\n   ERR invalid sort-permutation'); end
        sflag=1; sperm={sperm};
     end
     if sflag, A=sort(A,sperm{:}); end
  end

  sA=size(A); nA=prod(sA);
  i=find(sA>1); rA=numel(i);
  if rA<2
     if  sA(1)==1, sA(1)=[]; else rA=2; end
  end

  if isempty(nm)
     nm=inputname(1);
  end
  eflag=0;

  ov2s={'-f','nofac','sep'};

  cflag=1; nl=''; nl_='';
  if nA<=2 || vflag
     cflag=0; nl=char(10); if vflag>1, nl_='\n'; end
  end

  if isequal(nm,'ans'), nm=''; end

  if rA<2
       fmt=sprintf(' %%%dg',floor(log10(max(sA)))+1);
  else fmt=sprintf(',%%%dg',floor(log10(sA))+1); fmt(1)=' ';
  end
  if ~isempty(nm) && ~cflag, fmt(1)=[]; end
  fmt=regexprep(fmt,'%1g','%g');

  if nA==1
     if ~isempty(nm)
        if any(nm~=' ')
             s={sprintf('%s = ',nm)};
        else s={' '}; end
     else
        s={''};
     end
     display_1(A,m,Eflag,vflag,s{:});
  elseif nA>1
     ise=zeros(1,nA); nl2=nl;
     for i=1:nA, ise(i)=is_empty(A(i)); end

     for i=1:nA
        if i>1 && i<nA && all(ise(i-1:i+1))
           if ~eflag, fprintf(1,[nl_ '    :\n']); end % '\n...' : ┋┊
           eflag=eflag+1; continue;
        elseif eflag, eflag=0; if vflag<2, nl2=''; end
        else nl2=nl; end

        iA=ind2sub_aux(sA,i); l=sum(nm=='%');
        if l==length(iA), s=sprintf(nm,iA);
        elseif l==1,      s=sprintf(nm,i );
        else
           s=sprintf(fmt,iA);
           if cflag || isempty(nm), s=[s '. '];
           else s=[nm '(' s ') = ']; end
        end

        if ~ise(i)
           if ~cflag || vflag
                display_1(A(i),m,Eflag,vflag,s);
           else info(A(i),s,'-C'); end
        else fprintf(1,[nl2 '%s(empty)\n'],s); end
     end
  else
     s=sprintf('x%g',size(A));
     s={'empty QSpace array', s(2:end)};
     if ~isempty(nm)
          fprintf(1,[nl '   %s is %s (%s)\n' nl],nm,s{:});
     else fprintf(1,[nl '   (%s; %s)\n'      nl],   s{:});
     end
  end

  if ~cflag && nA
     q=[isempty(A(end).Q) isempty(A(end).data)];
     if ~q(1) || all(q), fprintf(1,'\n'); end
  end

end

% -------------------------------------------------------------------- %

function i=is_empty(A)
   if isempty(A) || isempty(A.data), i=1; else i=0; end
end

% -------------------------------------------------------------------- %
function kk = ind2sub_aux(s,k)

   if numel(s)<=2 && s(1)<2, kk=k; return; end

   k=k-1; kk=zeros(size(s));

   for i=1:length(s)
   kk(i)=mod(k,s(i)); k=(k-kk(i))/s(i); end

   kk=kk+1;

end

% -------------------------------------------------------------------- %
function display_1(A,m,Eflag,vflag,varargin)

  q=str2num(getenv('WB_VERBOSE'));
  if ~isempty(q) && q>=8
     wbrat(1); q=QSpace; normQS(q); getDimQS(q);
  end

  info(A,varargin{:});
  r=length(A.Q); if ~r, return, end

  dlen=numel(A.data);
  if ~dlen && ~isempty(A.Q) && isempty(A.Q{1})
     fprintf(1,'   (empty data and Q{:})\n');
     return;
  end

  cgflag=gotCGS(A); 
  sout=cell(1,dlen); l=1;

  if cgflag
     ns=length(find(A.info.qtype==','))+1;
     if ~isempty(A.Q) && ~isequal(size(A.info.cgr),[dlen, ns])
     error('Wb:ERR','CG size mismatch'); end
     nq=size(A.info.cgr,2);
     sfmt=sprintf('%%-%ds', max(10,3*length(A.Q)));
  else
     sfmt=sprintf('%%-%ds', 5+4*length(A.Q));
  end

  fstr={ '%11.4g'
         '  %-11s  %s%s\n'
         '  %s\n' }; % '  %8s\n'

  lfmt=isequal(get(0,'Format'),'long');
  if lfmt, fstr{1}='%16.13g'; end

  if ~isempty(A.Q)
     QQ=permute(cat(3,A.Q{:}),[3 2 1]);
     if ~isempty(A.info), t=A.info.qtype; else t=''; end
     qfmt=get_q_fmt(t,numel(A.Q),size(A.Q{1},2));
  else QQ=[]; end

  ov2s={'-f','nofac','sep'};

  m(end+1:2)=m(1); m(end+1:3)=0;
  mx=sum(m); mx=max(ceil(1.2*mx), mx+4); nx=0; xflag=0;

  ss=ones(dlen,r); sp=ones(dlen,1);
  for i=1:dlen
     q=size(A.data{i});
     ss(i,1:numel(q))=q; sp(i)=prod(q);
  end
  if r==1, ss=sort(ss,2);
     if all(ss(:,1)==1), ss=ss(:,2:end);
     else error('Wb:ERR','invalid rank-1 QSpace'); end
  end
  [~,is]=sort(sp);

  if dlen>max(4,4*mx) && m(3)
       ix=sort(is(end-m(3)+1:end));
  else ix=[]; end

  i=0;
  while i<dlen, i=i+1;
      if i>=m(1) && dlen>mx
         if ~isempty(ix), ix(find(ix<i))=[];
            if ~isempty(ix)
               if ix(1)-i<3
                    m(1)=ix(1);
               else i=ix(1); nx=nx+1; if xflag<1, xflag=1; end
               end
            end
         end
         if i>=m(1) && isempty(ix), j=dlen-m(2)+1;
            if j>i+2 && i<=dlen
               i=j; if nx~=1, xflag=1; end
            end
         elseif ~isempty(ix), ix(1)=[];
         end
         if xflag==1, sout{l}='     :   ...\n'; l=l+1; xflag=2; end
      end

      Ai=A.data{i}; 
      if r~=1, sa=size(Ai);
      else sa=ss(i,:); end

      s1=dim_to_str(sa,r);

      if cgflag, sc=cell(1,nq);
         for j=1:nq, sc{j}=cgr_size(A,i,j); end
         sc=cat2(1,sc{:},{1}); sc(:,end+1:r)=1; sa(end+1:r)=1;

         if ~vflag
            sc=prod(sc,1);
            s2=sprintf(sfmt,dim_to_str(sc,r));
         else
            n=size(sc,1); s2=cell(1,n);
            for j=1:n
               s2{j}=vec2str(sc(j,:),'fmt','%g',ov2s{:},'x');
            end
            s2=sprintf(' %6s',s2{:}); s2=s2(2:end);
         end

         sout{l}=sprintf(['%6d.  ' sfmt ' | %s' ],i,s1,s2);
      else
         sout{l}=sprintf(sfmt, sprintf('%6d.  %s',i,s1));
      end

      if ~isempty(QQ)
           sout{l+1}=[' [ ' sprintf(qfmt,QQ(:,:,i)') ' ]']; l=l+1;
      else sout{l+1}= ' [ ]'; end
      l=l+2;

      s=prod(sa); dfac=1;
      if 1 || Eflag
         if cgflag
            [dfac,sc]=get_cgr_fac(A,i,'-S');
            if ~isempty(sc); sc=['{' sc '}']; end
         else sc=''; end
      end

      if s==1, Ai=dfac*Ai;
         if isreal(Ai)
            q=sprintf(fstr{1}, Ai);
            if isempty(find(q=='.' | q=='e',1))
               q=[q '.']; if q(1)==' ', q(1)=[]; end
            end
            str=sprintf('  %s',q); if ~isempty(sc), str(end+1)=' '; end
         else
            str=num2str2(Ai,'fmt',fstr{1});
            str=sprintf('  %s',str);
         end
         if isempty(sc)
              sout{l}=[str '\n'];
         else sout{l}=[str ' ' sc '\n'];
         end; l=l+1;

         continue

      elseif Eflag
         if ~isreal(Ai), error('Wb:ERR',...
           '\n   ERR real data{} expected (got complex)'); end
         if Eflag==1, Ai=diag(Ai); sa=size(Ai); end
         if numel(find(sa>1))>1, error('Wb:ERR',...
           '\n   ERR invalid data{} (diagonal representation expected)'); end
         n=numel(Ai);
         if n<4
            q=sprintf([' ' fstr{1}],Ai);
         else q=[ ...
            sprintf([' ' fstr{1}],Ai(1:2)), sprintf(' ..(%d).. ',n-3), ...
            sprintf(fstr{1},Ai(end)) ];
         end
         sout{l}=sprintf(';  %s %s\n\n', q(2:end), sc); l=l+1;
         continue
      end

      if s, q=Ai(1); q=whos('q'); s=s*q.bytes;
      end

      if     s<2^10, s=sprintf('%6g b ',s);
      elseif s<2^20, s=sprintf('%6.1f k',s/2^10);
      else           s=sprintf('%6.1f M',s/2^20); end

      if ~isempty(sc)
           sout{l}=sprintf(fstr{2},s,sc,'');
      else sout{l}=sprintf(fstr{3},s);
      end; l=l+1;
  end

  fprintf(1,[sout{:}]);

end

% -------------------------------------------------------------------- %
function s=dim_to_str(sz,r)

   l=length(sz);
   if l<=r
      s=sprintf('x%g',sz); s=s(2:end);
   elseif l==r+1
      s=sprintf('x%g',sz(1:r));
      s=[ s(2:end) sprintf(' @%g',sz(r+1)) ];
   else
      s={ sprintf('x%g',sz(1:r)), sprintf('x%g',sz(r+1:end)) };
      s=[ s{1}(2:end), ' @' s{2}(2:end) ];
   end
end

% -------------------------------------------------------------------- %
