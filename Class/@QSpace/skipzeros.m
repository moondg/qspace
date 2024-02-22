function [A,n]=skipzeros(A,varargin)
% function A=skipzeros(A [,eps,'--all'])
%
%    Skip records in A with negligible weight (based on Frobenius norm).
%    See skipZerosQS for more information.
%
% Wb,Nov10,07

  n=zeros(size(A));

  for k=1:numel(A)
     [Ak,n(k)]=skipZerosQS(A(k),varargin{:});

     A(k)=QSpace(Ak);
  end

end

