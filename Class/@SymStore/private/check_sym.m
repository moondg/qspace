function [sym,r]=check_sym(sym)
% function [sym,r]=check_sym(sym)
% Wb,Aug01,19

  if isstruct(sym) && isfield(sym,'type')
     sym=sym.type;
  end

  if numel(sym)>2
     if isempty(regexp(sym,'^S[UOp]\d+$'))
        sym=regexprep(sym,'[()]','');
        if isempty(regexp(sym,'^S[UOp]\d+$'))
           wbdie('ERR invalid symmetry ''%s''',sym);
        end
     end
  elseif ~isequal(sym,'A4')
     wbdie('ERR invalid symmetry ''%s''',sym);
  end

  if nargout>1, r=get_rank(sym); end

end

