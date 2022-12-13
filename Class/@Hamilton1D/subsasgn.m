function HAM=subsasgn(HAM,S,val)
% function HAM=subsasgn(HAM,S,val)
% Reference http://www.cs.ubc.ca/~murphyk/Software/matlabObjects.html
% Classes and objects in matlab: the quick and dirty way
% Kevin Murphy, 2005

  HAM = builtin('subsasgn', HAM, S, val);

end

