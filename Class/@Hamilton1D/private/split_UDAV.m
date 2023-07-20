function UD=split_UDAV(AK,HD)
% function UD=split_UDAV(AK,HD)
% Wb,Aug11,16

% outsourced from update_psi_2site.m // Wb,Dec20,21

  if numel(AK.Q)~=2 || ~isequal(AK.Q{1},AK.Q{2})
     wbdie('invalid usage (AK)'); end
  if numel(HD.Q)~=2 || ~isequal(HD.Q{1},HD.Q{2})
     wbdie('invalid usage (HD)'); end

  [i1,i2,I]=matchIndex(AK.Q{1},HD.Q{1});
  if ~isempty(I.ix1), wbdie('invalid AK data !?'); end
  if ~isequal(i1,i2) || ~isempty(I.ix2)
     AK=struct(getsub(QSpace(AK),i1));
     HD=struct(getsub(QSpace(HD),i2));
  end

  nd=numel(HD.data);
  ss=zeros(nd,2); D1=cell(nd,1); Dc=cell(nd,1);
  for i=1:nd
     ss(i,:)=size(HD.data{i});

     n=ss(i,1); d1=zeros(n,1);
     for j=1:n, d1(j)=size(HD.data{i}{j,1},1); end
     D1{i}=d1; Dc{i}=cumsum(D1{i});
  end
  if norm(diff(ss,1,2)), wbdie('unexpected HD data !?'); end
  ndav=max(ss(:,1));

  d1=zeros(nd,1);
  for i=1:nd
     q=AK.data{i}; s=size(q); j=find(Dc{i}<=s(1));
     AK.data{i}=mat2cell(q,D1{i}(j),s(2));
     d1(i)=numel(j);
  end

  UD=QSpace(ndav,1);

  A_=QSpace(AK);

  for k=1:ndav
     X=struct(getsub(A_,find(d1>=k)));
     for j=1:numel(X.data), X.data{j}=X.data{j}{k,1}; end
     UD(k)=QSpace(X);

     t=UD(k).info.itags;
       t{1}=regexprep(t{1},'^[\d\w]+',sprintf('dav%g',k-1));
       t{2}=regexprep(t{2},'K\*','*');
     UD(k).info.itags=t;
  end

end

