function s=get_mstr(A)
% function s=get_mstr(A)
% outsourced from mtimes.m
% Wb,Dec05,11

   if ~isa(A,'SymOp') || numel(A)~=1
   error('Wb:ERR','\n   ERR invalid usage'); end

   s=A.istr;
   s=regexprep(s,'^[^=]+=[ ]+','');

   i=regexp(s,'[ .*+\-()\[\]]');
   if ~isempty(i) && (s(1)~=('(') || s(end)~=(')')), s=['(' s ')']; end

end

