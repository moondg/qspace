function M=full(A)
% wrapper routine for mpsFullQS()

  M=mpsFullQS(A(1));

  for i=2:length(A)
  M(:,:,i)=mpsFullQS(A(i)); end

return

