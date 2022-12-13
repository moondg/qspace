function M=full2(A)
% wrapper routine for mpsFull2QS()

  M=mpsFull2QS(A(1));

  for i=2:length(A)
  M(:,:,i)=mpsFull2QS(A(i)); end

return

