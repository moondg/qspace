function i=ne(A,B)
% function i=ne(A,B)
%
%    Simple wrapper to QSpace/isequal.m (whether two QSpaces are equal)
%    that overloads Matlab's the ~= (not-equal) operator.
%
% Wb,Feb22,24

  i=~isequal(A,B);

end

