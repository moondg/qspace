function [sym,r]=check_sym(sym)
% function [sym,r]=check_sym(sym)
% Wb,Aug01,19

  if isstruct(sym) && isfield(sym,'type')
     sym=sym.type;
  end

  if isempty(regexp(sym,'^S[UpO]\d+$'))
     sym=regexprep(sym,'[()]','');
     if isempty(regexp(sym,'^S[pON]\d+$'))
        wberr('ERR invalid symmetry ''%s''',sym);
     end
  end

  if nargout>1, r=get_rank(sym); end

end

