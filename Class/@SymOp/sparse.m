function A=sparse(A,eps)
% function B=sparse(A,eps)
% Wb,Dec05,11

  if nargin<2, eps=0; end
  for i=1:numel(A)
     if eps
        A(i).op(find(abs(A(i).op)<eps))=0;
        A(i).hc(find(abs(A(i).hc)<eps))=0;
     end

     if ~issparse(A(i).op), A(i).op=sparse(A(i).op); end
     if ~isempty(A(i).hc) && ~issparse(A(i).hc)
        A(i).hc=sparse(A(i).hc);
     end
  end

end

