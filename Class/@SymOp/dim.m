function s=dim(A)
% function s=dim(A)
% Wb,Dec05,11

   if ~isa(A,'SymOp'), error('Wb:ERR','\n   ERR invalid usage'); end

   s=size(A(1).op);
   if numel(s)>3 || s(1)~=s(2)
      error('Wb:ERR','\n   ERR invalid symmetry op'); end

   for i=2:numel(A)
       if ~isequal(s,size(A(i).op))
       error('Wb:ERR','\n   ERR severe size inconsistency'); end
   end

   s=s(1);

end

