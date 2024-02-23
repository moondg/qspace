function r=rank(A)
% function r=rank(A)
%
%    Returns the tensor-rank of QSpace A, i.e., its number of legs.
%
%  Wb,Mar02,08

  r=zeros(size(A));
  for i=1:numel(A), r(i)=length(A(i).Q); end

end

