function wbstack(k)
% Function wbstack
% Wb,Jan30,08

  stack=dbstack;

  if nargin, k=k+2; else k=2; end
  dispstack(stack(k:end));

end

