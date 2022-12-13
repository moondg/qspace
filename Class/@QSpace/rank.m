function r=rank(A)
% overloading rank function

  r=zeros(size(A));
  for i=1:numel(A), r(i)=length(A(i).Q); end

end

