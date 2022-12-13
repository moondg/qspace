function i=get_stype(HAM,k)
% function i=get_stype(HAM,k)
% Wb,Dec01,20

  if isfield(HAM.mpo,'stype'), i=HAM.mpo(k).stype;
  elseif isfield(HAM.info,'stype'), i=HAM.info.stype; else i=[]; end

  n=numel(i); if n~=1
     if ~n, i=1;
     else n(2)=numel(HAM.mpo);
        if diff(n), wberr('invalid HAM.info.stype (len=%g/%g)',n); end
        i=i(k);
     end
  end

end

