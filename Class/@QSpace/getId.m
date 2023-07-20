function A=getId(A,qs,k)
% function A=getId(A,qs [,k])
%
%    reduce QSpace A to rank-2 Id tensor which only has
%    vacuum state / scalar representation with each index.
%    if k is specified, itag{k} will be used.
%
% Wb,Oct18,14

% adapated from QSpace/getvac.m // Wb,Oct02,15

  if isempty(A),  wbdie('got empty QSpace'); end
  if numel(A)~=1, wbdie('single QSpace required'); end

  A=getsub(A,1); r0=numel(A.Q);
  A.data{1}=1;

  qs=reshape(qs,1,[]); nq=numel(qs); l0=0;

  if     r0>=2, A.Q=A.Q(1:2);
  elseif r0==1, A.Q=A.Q([1 1]);
  else wbdie('invalid usage'); end

  for i=1:2, A.Q{i}(:)=0; end

  for i=1:size(A.info.cgr,2), q=A.info.cgr(i);
     if ~isempty(q.type) && ~isempty(q.qset) && ~isempty(q.qdir)
       c=emptystruct(q);
       c.type=q.type;

       l=length(q.qset)/r0; l0=l0+l;
       c.qset=zeros(1,2*l);
       if l0<=nq, j=l0-l+1:l0;
          c.qset=[qs(j), qs(j)];
          A.Q{1}(1,j)=qs(j);
          A.Q{2}(1,j)=qs(j);
       end

       c.qdir='+-';
       c.cgw=['1' 0];

     A.info.cgr(i)=c; end
  end

  if nq && l0~=nq
     wbdie('invalid usage (length mismatch of qs: %d/%d)',nq,l0);
  end

  A.info.otype='';

  t=regexprep(A.info.itags,'*','');
  if nargin<3, k=0;
  elseif k>numel(t)
     error('Wb:ERR','\n   ERR index out of bounds (k=%d/%d)',k,numel(t));
  end

  if k, A.info.itags={t{k}, [t{k} '*']};
  else  A.info.itags={'', '*'}; end

end

