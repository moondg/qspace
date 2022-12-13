function Ek=get_ops_E(HAM,k,tflag)
% function Ek=get_ops_E(HAM,k,tflag)
% Wb,Apr22,14 ; Wb,Aug26,15

  if nargin<2 || nargin>3 || nargin==3 && ~isequal(tflag,'-t')
     helpthis, if nargin || nargout, wberr('invalid usage'), end
     return
  end
  if nargin>2
     if isequal(tflag,0) || isequal(tflag,'0'), tflag=0;
     else tflag=1; end
  end

  i=get_stype(HAM,k);
  Ek=HAM.oez(1,i).op;

  if tflag
     if ~isIdentityQS(Ek), wberr(...
        'invalid HAM.mpo(%g): expecting HAM.oez(1,%g) to be Id',k,i);
     end
     for i=1:numel(Ek.data)
        q=Ek.data{i}; rq=round(q); e=norm(q-rq);
        if e && e<1E-12, Ek.data{i}=rq; end
     end
  end

end

