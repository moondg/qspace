function A=full(A)
% function B=full(A)
% Wb,Dec05,11

  for i=1:numel(A)
     if issparse(A(i).op), A(i).op=full(A(i).op); end
     if ~isempty(A(i).hc) && issparse(A(i).hc)
        A(i).hc=full(A(i).hc);
     end
  end

end

