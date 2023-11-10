function s=dim(A)
% function s=dim(A)
% Wb,Dec05,11

   if ~isa(A,'SymOp'), wbdie('invalid usage'); end

   s=size(A(1).op);
   if numel(s)>3 || s(1)~=s(2)
      wbdie('invalid symmetry op'); end

   for i=2:numel(A)
       if ~isequal(s,size(A(i).op))
       wbdie('severe size inconsistency'); end
   end

   s=s(1);

end

