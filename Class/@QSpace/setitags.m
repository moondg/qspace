function varargout=setitags(A,varargin)
% function A=setitags(A,...)
%
%   set itag(s) while preserving conj flag.
%
% Usage #1: A=setitags(A,{t1,t2,...} [,[i1,i2,...]])
%
%   set specified or all itags
%   setitags(A,{'K1','K2','s2'})        sets itags { 'K1','K2,'s2' }
%   setitags(A,{'K','K','s'},{1,2,2})   same (numbers are appended to itags!)
%
% Usage #2: A=setitags(A,tt,i)
%
%   set specific itags tt to itag{i}
%
% Usage #3: set itags for A-tensor (conj-flags always remain in place)
%
%   A=setitags(A,'-A',k)           {'K<k-1>','K<k>','s<k>}
%   A=setitags(A,'-A:K,K,s',  k)   {'K<k-1>','K<k>','s<k>} more specific
%   A=setitags(A,'-A>K,K,s@n',k)   same, also using n-digit format (default: 2)
%   A=setitags(A,'-A<K,K,s@n',k)   {'K<k>','K<k+1>','s<k>} i.e. R->L order
%
% Usage #4: set itags for operator
%
%   S=setitags(S,'-op:X', k )       {'<X><k>','<X><k>','op'}
%   S=setitags(S,'-op:X','s')       {'<X>','<X>','s'}
%   S=setitags(S,'-op:<regex>', A)  set itags based on behavior
%                                   that mimicks contractQS
% Usage #5
%
%   A=setitags(A,ia,B,ib);  copy itags from B to A
%
% Usage #6: any of the usages above, but return itags rather than A,
% by requesting 2 output args, e.g.
%
%   [tt,i]=setitags(A,...)
%   [tt,~]=setitags(A,...)
%
%   this returns itags+indizes themselves, rather than setting the
%   itags in A (the conjugation is inherited from A; hence only a
%   single QSpace A is expected).
%
% Wb,Mar24,16 ; Wb,May28,18

  if nargin<2
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  idx=[]; n=numel(varargin{1});

  if nargin==4 && (isstruct(varargin{2}) || isa(varargin{2},'QSpace')) ...
     && isnumeric(varargin{1}) && isnumeric(varargin{3})

     ia=varargin{1}; ib=varargin{3};n=numel(ia);
     if numel(ib)~=n
        wberr('invalid usage (length mismatch (%g/%g)',n,numel(ib));
     end

     tB=varargin{2}.info.itags; tB=regexprep(tB,'\**$','');

     if ~n, n=[numel(tA), numel(tB)];
        if diff(n)
           wberr('invalid usage (length mistmatch %g/%g)',n); end
        n=n(1); ia=1:n; ib=1:n;
     end

     idx=ia; tt=tB(ib); w=5;

  elseif ischar(varargin{1}), q=varargin{1};
     if ~isempty(regexp(q,'^-A'))
        if nargin~=3 || ~isnumber(varargin{2}), error('Wb:ERR',...
           '\n   ERR invalid usage (missing or invalid site index k)');
        end
        k=varargin{2}; fmt='%02g'; q_=q;

        i=regexp(q,'@\d+$');
        if ~isempty(i)
           fmt=sprintf('%%0%gg',str2num(q(i+1:end)));
           q=q(1:i-1);
        end

        if numel(q)<=2, t={'K','K','s'; k-1, k, k};
        else
           if ~isempty(regexp(q,'^-A[:<>][a-zA-Z,]+$'))
                t=textscan(q(4:end),'%s','whitespace',' ,');
           else t={''}; end
           if numel(t)~=1 || numel(t{1})~=3, error('Wb:ERR',...
             '\n   ERR invalid usage (q=''%s'' %g/3 !?)',q_,numel(t{1}));
           end

           t=reshape(t{1},1,[]);
           if     q(3)=='<', t(2,:)={k+1 k, k};
           elseif q(3)=='>', t(2,:)={k k+1, k};
           else              t(2,:)={k-1 k, k};
           end
        end

        fmt=['%s' fmt];
        for i=1:3, t{1,i}=sprintf(fmt,t{1,i},t{2,i}); end
        tt=t(1,:);

        idx=1:3; n=3; w=3;

     elseif n>=3 && isequal(q(1:3),'-op')
        if n>3 && q(4)==':', q=q(5:end);
        else wberr('invalid usage'); end

        i=''; q=regexprep(q,'@(\d+)$(?@i=$1;)','');
        if ~isempty(i)
             fmt=sprintf('%%0%gg',str2num(i));
        else fmt='%02g'; end

        if nargin>2, k=varargin{2}; else k=[]; end
        if isa(k,'QSpace') || (isfield(k,'Q') && isfield(k,'info'))
           k=k.info.itags; q=split(q,':');
           if ~isempty(regexp(q{1},'[\^\$\[\]\*\+]'));
               i=findstrc(k,q{1}); n=numel(i);
               if n==1, tt={ k{i}, k{i}, ''};
               elseif ~n, wberr('failed to find matching itag'); 
               else wberr('multiple matching itags found'); end
           else
               tt={ q{1}, q{1}, ''};
           end
           n=numel(q); if n>1, tt(3:n+1)=q(2:end); end
        else
           if ~isempty(q), tt={q,q,'op'}; else tt={'','',''}; end
           if nargin>2, fmt=['%s' fmt];
              if isnumber(k), for i=1:2, tt{i}=sprintf(fmt,tt{i},k); end
              elseif ischar(k), tt{3}=k;
              else wberr('invalid usage'); end
           end
        end
        idx=1:3; n=3; w=4;
     end
  end

  if isempty(idx)
     tt=varargin{1}; if ~iscell(tt), tt={tt}; end
     if nargin>2
        idx=varargin{2}; n=[ numel(idx), numel(tt)];
        if diff(n), display(idx), wberr('invalid indices (%g/%g)',n);
        else n=n(1); end
     else
        n=numel(tt);
     end

     for i=1:n
        if ~ischar(tt{i}), display(tt{i}), wberr('invalid itag'); end
        tt{i}=regexprep(tt{i},'\**$','');
        if ~isempty(regexp(tt{i},'[\s\*,;()\-]')) || length(tt{i})>8
           display(tt{i}), wberr('invalid itag'); end
     end

     if ~isempty(idx) && isnumeric(idx), w=2;
     else w=1;
        if ~isempty(idx)
           for i=1:n, tt{i}=sprintf('%s%g',tt{i},idx{i}); end
        end
        idx=1:n;
     end
  end

  nidx=numel(idx);
  if nidx, q=unique(idx); q=[ numel(q), nidx, q([1 end]) ];
     if diff(q(1:2)), idx, wberr('invalid usage (non-unique index)'); end
     if q(3)<1 || q(4)>32, idx, wberr('(likely) invalid index'); end
  else wblog('WRN','got empty idx (%g/%g)',nidx,numel(tt));
  end

  nA=numel(A); TK=cell(1,nA);
  for k=1:nA
     if ~isfield(A(k).info,'itags'), continue; end
     tk=A(k).info.itags;

     nk=numel(tk);
     if w==4, if nk<2 || nk>4, wberr(...
       'unexpected rank-%g operator',nk); end
     elseif any(idx>nk), wberr(...
       'invalid usage (index out of bounds %g/%g)',max(idx),nk); 
     end

     for i=1:nidx, j=idx(i);
        if j<=nk % skip 'op' if scalar operator // Wb,Jul24,22
           tk{j}=set_itag(tk{j},tt{i});
        end
     end
     TK{k}=tk;
  end

  if nargout<=1
     for k=1:nA, A(k).info.itags=TK{k}; end
     if ~nargout, n=inputname(1);
        if ~isempty(n), assignin('caller',n,A); return; end
     end
     varargout={A};
  else
     if numel(A)~=1, error('Wb:ERR','\n   ERR invalid usage'); end
     varargout={TK{1},idx};
  end

end

% -------------------------------------------------------------------- %

function t=set_itag(t,tref)

   t=regexprep(t,'^[^\*]*',tref,'emptymatch');

end

% -------------------------------------------------------------------- %

