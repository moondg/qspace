function [s]=issingleop(A)
% function [s]=issingleop(A [,opts])
%
%    check for possibly non-scalar, yet still rank-2 operators.
%
% See also @QSpace/isscalarop for scalar, i.e. block-diagonal
% (rank-2) operators.
%
% Wb,May11,13

  s=ones(size(A));
  for k=1:numel(A), r=numel(A(k).Q);
     if r<2 || r>3 || r==3 && dim(A(k),3,'-f')~=1, s(k)=0; end
  end

end

