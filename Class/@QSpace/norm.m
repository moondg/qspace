function M=norm(A)
% function M=norm(A)
%
%    Simple wrapper routine for normQS().
%
% Wb,Oct18,14

  M=zeros(size(A));
  for i=1:numel(A), M(i)=normQS(A(i)); end

end

