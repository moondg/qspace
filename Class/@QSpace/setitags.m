function varargout=setitags(A,varargin)
% function A=setitags(A,...)
%
%   set itag(s) while preserving conj flag with usage #1 - #6
%   as detailed below. Note that conjugate flags always remain
%   preserved, i.e., are ignore in the input arguments if present.
%
% Usage #1: A=setitags(A, [{i1,i2,...},] tt)    ^1)
%
%   set specified leading set of, or all itags tt = {t1,t2,...}
%   setitags(A,{'K','K','s'},{1,2,2})   using cell, numbers are appended to itags
%   setitags(A,{'K1','K2','s2'})        same as above, specifying explicit itags
%
% Usage #2: A=setitags(A,i,tt)    ^1)
%
%   set itags for indices i to specified itags tt = {t1,t2,...}
%   (note that this differentiates from usage #1 by using a
%   numeric index array i, rather than a cell array with usage #1).
%
% Usage #3: set itags for A-tensors
%
%   MPS-like basis transformations with index order convention
%   LRs = (left,right, local state space) in mind
%
%   A=setitags(A,'-A',k)           {'K<k-1>','K<k>','s<k>}
%   A=setitags(A,'-A:K,K,s',  k)   {'K<k-1>','K<k>','s<k>} more specific
%   A=setitags(A,'-A>K,K,s@n',k)   same, also using n-digit format (default: 2)
%   A=setitags(A,'-A<K,K,s@n',k)   {'K<k>','K<k+1>','s<k>} i.e. R->L order
%
% Usage #4: set itags for operators
%
%   with index order convention X,X',op = [bra,ket,irop index]
%
%   S=setitags(S,'-op:X','op')       {'<X>','<X>','op'}
%   S=setitags(S,'-op:X[,op[#]][@n]',k)    {'<X><k>','<X><k>','op[k]'}
%     - where the default operator itag if not specified is 'op'
%     - if operator itag contains trailing #
%       the index k will be added in the same format as with X
%     - the optional trailing '@n' specifies n-digit format (default: n=2)
%
%   S=setitags(S,'-op:<regex>', A)  set itags based on behavior
%                                   that mimicks contractQS
% Usage #5
%
%   A=setitags(A,ia,B,ib);  copy itags from B to A
%
% Usage #6: any of the usages above, but return the final itags
% rather than A, by requesting 2 output args, e.g.
%
%   [tt,i]=setitags(A,...)
%   [tt,~]=setitags(A,...)
%
%   this returns itags+indizes themselves, rather than setting the
%   itags in A (the conjugation is inherited from A; hence only a
%   single QSpace A is expected).
%
% ^1) the order of arguments for usage #1 and #2 changed
%   switching the 2nd and 3rd argument for more intuitive reasons;
%   therefore for backward compatibility [07/05/2023]:
%   usage #1 and usage #2 also support the reverse order of arguments
%   2 and 3, i.e., permits one to specify the index (set) i last.
%
% Wb,Mar24,16 ; Wb,May28,18

  if nargin<2
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  idx=[]; n=numel(varargin{1});

  if nargin==4 && (isstruct(varargin{2}) || isa(varargin{2},'QSpace')) ...
     && isnumeric(varargin{1}) && isnumeric(varargin{3})
     w=5;

     ia=varargin{1}; ib=varargin{3}; n=numel(ia);
     if numel(ib)~=n, wbdie(...
       'invalid usage (index length mismatch %g/%g)',n,numel(ib));
     end

     tB=varargin{2}.info.itags; tB=regexprep(tB,'\**$','');

     if ~n, n=[numel(tA), numel(tB)];
        if diff(n)
           wbdie('invalid usage (length mistmatch %g/%g)',n); end
        n=n(1); ia=1:n; ib=1:n;
     end

     idx=ia; tt=tB(ib);

  elseif ischar(varargin{1}), q=varargin{1};
     if ~isempty(regexp(q,'^-A'))
        midx=[]; w=3;

        if nargin<3 || ~isnumber(varargin{2})
           wbdie('invalid usage (missing or invalid site index k)');
        elseif nargin>3
           getopt('init',varargin(3:end));
              midx=getopt('--mark',[]);
           getopt('check_error');
        end
        k=varargin{2}; fmt='%02g'; q_=q;

        i=regexp(q,'@\d+$');
        if ~isempty(i), n=str2num(q(i+1:end));
           if n~=2, fmt=sprintf('%%0%gg',n); end
           q=q(1:i-1);
        end

        if numel(q)<=2, t={'K','K','s'; k-1, k, k};
        else
           if ~isempty(regexp(q,'^-A[:<>][a-zA-Z,]+$'))
                t=textscan(q(4:end),'%s','whitespace',' ,');
           else t={''}; end
           if numel(t)~=1 || numel(t{1})~=3
              wbdie('invalid usage (q=''%s'' %g/3 !?)',q_,numel(t{1}));
           end

           t=reshape(t{1},1,[]);
           if     q(3)=='<', t(2,:)={k+1 k, k};
           elseif q(3)=='>', t(2,:)={k k+1, k};
           else              t(2,:)={k-1 k, k};
           end
        end

        fmt=['%s' fmt];
        for i=1:3, t{1,i}=sprintf(fmt,t{1,i},t{2,i}); end

        if ~isempty(midx)
            if any(idx>3) || any(diff(sort(idx)<1))
               wbdie('invalid usage (invalid --mark index)'); end
            for i=midx, t{1,i}=[ t{1,i} '''' ]; end
        end
        tt=t(1,:);

        idx=1:3; n=3;

     elseif ~isempty(regexp(q,'^-op'))
        w=4;

        if n<4 || q(4)~=':', wbdie('invalid usage'); end
        i=0; t=regexprep(q(5:end),'@(\d+)$(?@i=str2num($1);)','');
        if i && i~=2, fmt=sprintf('%%0%gg',i);
        else fmt='%02g'; end

        if nargin>2, k=varargin{2}; else k=[]; end
        if isa(k,'QSpace') || (isfield(k,'Q') && isfield(k,'info'))
           k=k.info.itags; t=split(t,':');
           if ~isempty(regexp(t{1},'[\^\$\[\]\*\+]'));
               i=findstrc(k,t{1}); n=numel(i);
               if n==1, tt={ k{i}, k{i}, ''};
               elseif ~n, wbdie('failed to find matching itag'); 
               else wbdie('multiple matching itags found'); end
           else
               tt={ t{1}, t{1}, ''};
           end
           n=numel(t); if n>1, tt(3:n+1)=t(2:end); end
        else
           kstr=ischar(k) && ~isempty(k);
           if ~isempty(t)
              t=textscan(t,'%s','whitespace',' ,'); t=t{1}; n=numel(t);
              if     n==1, tt={t{1},t{1},'op'};
              elseif n==2, tt=t([1 1 2]);
                 if kstr, wbdie(['invalid usage ' ...
                    '(multiple operator itags specified: %s / %s)'],q,k);
                 end
              else wbdie('invalid usage (got %d itags with -op)',n); end
           else tt={'','',''}; end

           if nargin>2
              if kstr, tt{3}=k;
              elseif isnumber(k), n=2; fmt=['%s' fmt];
                 if ~isempty(tt{3}) && tt{3}(end)=='#'
                    tt{3}=regexprep(tt{3},'#+$',''); n=3;
                 end
                 for i=1:n, tt{i}=sprintf(fmt,tt{i},k); end
              else wbdie('invalid usage'); end
           end
        end
        idx=1:3; n=3;
     end
  end

  if isempty(idx)
     tt=varargin{1}; if ischar(tt), w1=1; tt={tt};
     elseif iscell(tt) && ischar(tt{1}), w1=2; else w1=0; end

     if nargin>2
        if w1, idx=varargin{2}; 
        else
           idx=tt; tt=varargin{2};
           if ischar(tt), w1=-1; tt={tt}; end
        end

        n=[ numel(idx), numel(tt)];
        if diff(n), display(idx), wbdie('invalid indices (%g/%g)',n);
        else n=n(1); end
     else
        n=numel(tt);
     end

     for i=1:n
        if ~ischar(tt{i}), display(tt{i}), wbdie('invalid itag'); end
        tt{i}=regexprep(tt{i},'\**$','');
        if ~isempty(regexp(tt{i},'[\s\*,;()]')) || length(tt{i})>8
           wbdie('invalid itag (%s)',tt{i}); end
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
     if diff(q(1:2)), idx, wbdie('invalid usage (non-unique index)'); end
     if q(3)<1 || q(4)>32, idx, wbdie('(likely) invalid index'); end
  else wblog('WRN','got empty idx (%g/%g)',nidx,numel(tt));
  end

  nA=numel(A); TK=cell(1,nA);
  for k=1:nA
     if ~isfield(A(k).info,'itags'), continue; end
     tk=A(k).info.itags;

     nk=numel(tk);
     if w==4, if nk<2 || nk>4, wbdie(...
       'unexpected rank-%g operator',nk); end
     elseif any(idx>nk), wbdie(...
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
     if numel(A)~=1, wbdie('invalid usage'); end
     varargout={TK{1},idx};
  end

end

% -------------------------------------------------------------------- %

function t=set_itag(t,tref)

   t=regexprep(t,'^[^\*]*',tref,'emptymatch');

end

% -------------------------------------------------------------------- %

