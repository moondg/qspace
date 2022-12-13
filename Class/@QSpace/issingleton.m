function i=issingleton(A,kk)
% Function i=issingleton(A,kk)
% Wb,Apr24,08

  if nargin<2 || isequal(kk,'all'), kk=1:length(A.Q);
  else kk=reshape(kk,1,[]); end

  DD=datasize(A); i=zeros(size(kk));

  for j=1:length(kk), k=kk(j);
      if norm(diff(A.Q{k},1))~=0, continue; end
      if all(DD(:,k)==1), i(j)=1; end
  end

  if nargin<2, i=any(i~=0); end

end

