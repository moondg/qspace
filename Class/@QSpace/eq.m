function i=eq(A,B)
% function i=eq(A,B)
%
%    Simple wrapper to QSpace/isequal.m (whether two QSpaces are equal)
%    that overloads Matlab's the == (equal) operator.
%
% Wb,Feb22,24

  i=isequal(A,B);

end

