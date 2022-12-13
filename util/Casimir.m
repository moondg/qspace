function c2=Casimir(sym,q)
% function Casimir(sym,q)
% Wb,Oct05,15

  if nargin~=2 || size(q,1)~=1 || nargout>1
     error('Wb:ERR','\n   ERR invalid usage'); end

  if isnumeric(q), q=sprintf('%g',q); end

  c2=SymStore(sym,'Casimir',q);

end

