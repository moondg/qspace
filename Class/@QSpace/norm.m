function M=norm(A)
% wrapper routine for normQS()

  M=zeros(size(A));
  for i=1:numel(A), M(i)=normQS(A(i)); end

end

