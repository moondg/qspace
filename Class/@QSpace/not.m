function q=not(A)
% function q=not(A)
%
%    Overload of conversion to logicals (unary operator)
%    e.g. for conditionals: if ~A, ...; else ...; end
%
% Wb,Feb18,21

  q = ~logical(A);

end

