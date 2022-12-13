function s=length(P,vn)
% Function length(P,vn)
% get size by name (in order to keep PSet dynamic)
% Wb,Jul31,08

  if nargin<2
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

  for i=1:P.r
  if isequal(P.vars{i},vn), s=P.s(i); return; end, end

  error('Wb:ERR','var ''%s'' not contained within PSet',vn); 

end

