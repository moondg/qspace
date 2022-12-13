function C=kron(A,B)
% overloading kron()

  if isempty(A.Q) && isempty(A.data) C=QSpace(); return, end
  if isempty(B.Q) && isempty(B.data) C=QSpace(); return, end

  C=QSpace(mpsTensorProdQS(A,B));

return

