function c2=Casimir(S)
% function Casimir(S)
% Wb,Oct05,15

  if numel(S)~=1 || nargout>1
     error('Wb:ERR','\n   ERR invalid usage'); end

  [q,d,D]=getQDimQS(S,'op');

  if size(q,1)~=1, error('Wb:ERR',...
    '\n   ERR spin-operator in single local multiplet required'); end

  c2=SymStore(S.info.qtype,'Casimir',sprintf('%g',q));

end

