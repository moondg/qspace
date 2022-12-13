function A=subsasgn(A,S,val)
% SUBSASGN
%
% Reference http://www.cs.ubc.ca/~murphyk/Software/matlabObjects.html
% Classes and objects in matlab: the quick and dirty way
% Kevin Murphy, 2005

  if ~isa(A,  'QSpace'), A=QSpace(A); end

% required by matlab/2016a // Wb,Aug03,16
  if numel(S)==1 && isequal(S.type,'()') && ~isa(val,'QSpace')
     val=QSpace(val);
  end

  A = builtin('subsasgn',A,S,val);

end

