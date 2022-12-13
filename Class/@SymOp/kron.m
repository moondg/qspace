function C=kron(A,B)
% overloading kronecker product
% Wb,Nov30,11

  os=' (x) ';

  if isa(A,'SymOp') && isa(B,'SymOp'), C=A; n=numel(A);
     if n~=numel(B), error('Wb:ERR','\n   ERR size mismatch'); end
     for i=1:n
        C(i).istr=[ get_mstr(A(i)) os B(i).istr ];
        C(i).op=mykron(C(i).op,B(i).op);
        C(i).hc=mykron(C(i).hc,B(i).hc);
     end
  elseif isa(B,'SymOp') && isnumeric(A), C=B; s=num2str(A);
     for i=1:numel(C)
        C(i).istr=[ get_mstr(C(i)) os s ];
        C(i).op=mykron(A,C(i).op);
        C(i).hc=mykron(A,C(i).hc);
     end
  elseif isa(A,'SymOp') && isnumeric(B), C=A; s=num2str(B);
     for i=1:numel(C)
        C(i).istr=[ get_mstr(C(i)) os s ];
        C(i).op=mycron(B,C(i).op);
        C(i).hc=mycron(B,C(i).hc);
     end
  else error('Wb:ERR','\n   ERR invalid usage'); end

end

function AB=mykron(A,B)
    e=[isempty(A),isempty(B)];
    if any(e)
       if xor(e(1),e(2)), error('Wb:ERR',...
         '\n   ERR invalid usage (mixed empty / non-empty)');
       end
       AB=A; return
    end
    AB=kron(A,speye(size(B))) + kron(speye(size(A)),B);
end

