function Aout=subsref(A,S)
% SUBSREF
%
% Reference http://www.cs.ubc.ca/~murphyk/Software/matlabObjects.html
% Classes and objects in matlab: the quick and dirty way
% Kevin Murphy, 2005

  Aout=builtin('subsref',A,S);

end

