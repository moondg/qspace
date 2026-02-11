function C=contract(varargin);
% function C=contract(A,B [,ia,ib [,perm]]);
%
%   general routine to contract specified sets of indizes
%   ica and icb of the pair of tensors A and B, respectively.
%
% Usage 2: contract(A, [i1 i2; i3 i4; ...])
%
%   generalized trace
%   that traces out i1 with i2, i3 with i4, etc.
%
% Usage 3: redirect to contractQS()
%
%   in case of incompatible usage on bare numeric tensors
%   assume QSpace contracteion, instead, i.e., call contractQS().
%
% AWb, F.Verstraete

% [Wb,07/25/2023] added usage #3 (redirect to contractQS)

  nargs=nargin;
  if nargs
     for i=1:nargs, x=varargin{i};
        if isstruct(x) 
           if isfield(x,'Q') && isfield(x,'data')
              i=nargs+1; if numel(x)~=1, i=-i; end
              break
           end
        elseif iscell(x) && ~isnumeric(x)
           i=nargs+2; if numel(x)<2, i=-i; end
           break
        end
     end
     if abs(i)>nargs
        if i>0
           C=contractQS(varargin{:});
           C=QSpace(C); return
        end
        wbdie('invalid usage #3 (failed checks for contractQS)');
     end
  end

  if nargs<2, helpthis
     if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  A=varargin{1};
  B=varargin{2};

  if nargs==2, sa=size(A); ra=ndims(A); [n,m]=size(B);
     if isequal(B,round(B)) && m==2, i=sort(B(:)); 
        if ~n, wbdie('empty index set for generalized trace'); end
        if any(i<1 | i>ra) || any(diff(i)<1), wbdie(...
           'invalid index for generalized trace'); end
        if any(diff(sa(B),[],2)), wbdie(...
           'size mismatch for generalized trace'); end
        i2=sort(B,2);

        for i=1:n
           ia=i2(1,:);
           i1=ia(1); i2=ia(2);
           ka=1:ra; ka(ia)=[];
           sA=[prod(sa(ka)), prod(sa(ia))]; pA=[ka ia];
           A=reshape(permute(A,pA),sA)*reshape(eye(sa(i1)),[],1);

           if i<n
              ipA(pA)=1:ra;
              i2=ipA(i2(2:end,:));
           end

           sa=sa(ka); sa(end+1:2)=1; ra=length(sa);
           A=reshape(A,sa);
        end
        C=A;
     else
        if ~isequal(sa,size(B))
           wbdie('size mismatch for full contraction'), end
        C=reshape(A,1,[])*reshape(B,[],1);
     end
     return
  end

  if nargs<4 || nargs>5, wbdie('invalid usage'); end
  ia=varargin{3};
  ib=varargin{4};

  if iscell(ia), ra=ia{2}; ia=ia{1}; else ra=max(ia); end
  if iscell(ib), rb=ib{2}; ib=ib{1}; else rb=max(ib); end

  sa=size(A); sa(end+1:ra)=1; ra=length(sa);
  sb=size(B); sb(end+1:rb)=1; rb=length(sb);

  ka=1:ra; ka(ia)=[];
  kb=1:rb; kb(ib)=[];

  if isempty(ia)
     wbdie('invalid usage (empty contraction index set)'); 
  elseif ~isequal(sa(ia),sb(ib))
     wbdie('size mismatch of contracted indices: %s <> %s',...
     vec2str(sa(ia),'-f'), vec2str(sb(ib),'-f'));
  end

  sA=[ prod(sa(ka)), prod(sa(ia)) ]; pA=[ka ia];
  sB=[ prod(sb(ib)), prod(sb(kb)) ]; pB=[ib kb];
  sc=[sa(ka) sb(kb)]; sc(end+1:2)=1;

  A2=reshape(permute(A,pA),sA);
  B2=reshape(permute(B,pB),sB);
  C =reshape(A2*B2,sc);

  if nargs>4 && ~isempty(varargin{5}), pC=varargin{5};
     rc=numel(sc); l=numel(pC);
     if ~isperm(pC), wbdie('invalid input permutation for C'); end
     if l>rc, wbdie('invalid input permutation for C (length %d/%d)',l,rc); end
     if l<rc, pC(l+1:rc)=l+1:rc; end
     C=permute(C,pC);
  end

end

