function A=uminus(A)
% function A=uminus(A)
%
%    overloading the unary minus operator (as in -A)
%
% Wb,Sep17,06

  for k=1:prod(size(A))
     for i=1:length(A(k).data)
        A(k).data{i}=-A(k).data{i};
     end
  end

end

