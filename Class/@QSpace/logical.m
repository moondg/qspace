function q=logical(A)
% function q=logical(A)
%
%    Overload of conversion to logicals (unary operator)
%    e.g. for conditionals: if A, ...; else ...; end
%
% see also QSpace/not.m [which returns ~logical()]
% Wb,Feb18,21

  q=true(size(A));
  for i=1:numel(A)
      if isempty(A(i).Q) && isempty(A(i).data), q(i)=0; end
  end

end

