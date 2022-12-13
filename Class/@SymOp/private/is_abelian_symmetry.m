function i=is_abelian_symmetry(sym)
% function i=is_abelian_symmetry(sym)
% Wb,Aug29,16

  if isequal(sym,'A') || isequal(sym,'P') || ~isempty(regexp(sym,'^Z\d+$'))
       i=1;
  else i=0; end

end

