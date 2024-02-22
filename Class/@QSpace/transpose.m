function A=transpose(A)
% function A=transpose(A)
%
%   This function overloads the .' operator in Matlab
%   for `element wise transpose'.

% Wb,Aug19,06

  if isempty(A) || numel(A)==1 && isempty(A.Q) && isempty(A.data)
     return
  end

  A=builtin('transpose',A);

  for k=1:numel(A), r=numel(A(k).Q);
     if r==2, p=[2 1];
     elseif r==3 && isequal(A(k).info.otype,'operator'), p=[2 1 3];
     else wbdie('got rank-%d QSpace',r); end
     A(k)=QSpace(permuteQS(A(k),p));
  end

end

