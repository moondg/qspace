function C=acomm(A,B)
% function C=acomm(A,B)
% Wb,Nov17,11

  if isnumeric(A) && isa(B,'SymOp')

     C=B; if issparse(A), A=full(A); end
     for i=1:numel(C)
        if ~isempty(C(i).op)
             C(i).hc=acomm(A,C(i).op);
        else C(i).hc=[]; end
     end

  elseif isnumeric(B) && isa(A,'SymOp')

     C=A; if issparse(B), B=full(B); end
     for i=1:numel(C)
        if ~isempty(C(i).op)
             C(i).hc=acomm(C(i).op,B);
        else C(i).hc=[]; end
     end

  elseif isa(A,'SymOp') && isa(B,'SymOp')

     if numel(A)~=1 || numel(B)~=1, error('Wb:ERR',['\n   ' ...
    'ERR invalid usage (A: %s, B: %s)'],sizestr(A),sizestr(B)); end

     C=A; C.hc=[]; C.istr=['{' A.istr ', ' B.istr '}'];
     C.op=acomm(A.op,B.op);

  else error('Wb:ERR','\n   ERR invalid usage'); end

end

