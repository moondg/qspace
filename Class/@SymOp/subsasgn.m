function G=subsasgn(G,S,val)
% function G=subsasgn(G,S,val)
% Wb,Nov18,11

   if ~isa(G,'SymOp'), G=SymOp(G); end
   G=builtin('subsasgn',G,S,val);

end

