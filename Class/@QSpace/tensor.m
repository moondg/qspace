function C=tensor(varargin)
% function C=tensor(A,B,...)
% kron including merging of quantum numbers
% Wb,Jul07,09

  if nargin==1 && numel(varargin{1})>1
     n=numel(varargin{1}); ops=cell(1,n);
     for i=1:n, ops{i}=varargin{1}(i); end
  else ops=varargin; end

  nops=length(ops);
  r=repmat([-1 0],nops,1); sq=zeros(1,nops); fac=1;

  for i=1:nops, A=ops{i};
     if isempty(A) || isempty(A.Q)
        if ~isempty(A.data)
           if ~length(A.data)==1 || ~isscalar(A.data{1})
           error('Wb:ERR','invalid QSpace object'); end
           fac=fac*A.data{1};
        else C=QSpace; return, end
        continue
     end
     r(i,:)=[length(A.Q), size(A.Q{1},2)];
  end

  i=find(r(:,1)<0);
  if ~isempty(i) r(i,:)=[]; sq(i)=[]; ops(i)=[]; nops=length(ops); end

  if norm(diff(r,[],1)), r
  error('Wb:ERR','QSpace::tensor - rank or QSpace mismatch'); end
  if any(r(:,1)~=2)
  error('Wb:ERR','QSpace::tensor - require operators (rank-2)'); end

  C=ops{1};
  for i=2:nops, C=tensor_aux(C,ops{i}); end

  if fac~=1, for i=1:length(C.data)
  C.data{i}=fac*C.data{i}; end, end

end

% -------------------------------------------------------------------- %

function C=tensor_aux(A,B)

  C=QSpace(mpsTensorProdQS(A,B));

  ra=length(A.Q);
  rb=length(A.Q); if ra~=2 || rb~=2, error('Wb:ERR','invalid usage'); end

  qq={ expandQ(A,'op'), expandQ(B,'op') };
  E=QSpace(qq{:},'-Rlast','identity');

  C=QSpace(contractQS(E,1:2,contractQS(C,3:4,E,1:2),1:2));

end

% -------------------------------------------------------------------- %

function C=tensor_aux_not_quite(A,B)

  ra=length(A.Q);
  rb=length(B.Q); if ra~=2 || rb~=2, error('Wb:ERR','invalid usage'); end

  C=QSpace(mpsTensorProdQS(A,B));

  Q=C.Q; rc=length(Q); if rc~=4, error('Wb:ERR','???'); end
  Q{1}=Q{1}+Q{2}; Q{2}=Q{3}+Q{4};
  C.Q=Q(1:2);

  for i=1:length(C.data); d=C.data{i};
  s=size(d); C.data{i}=reshape(d,s(1)*s(2),[]); end

end

% -------------------------------------------------------------------- %

