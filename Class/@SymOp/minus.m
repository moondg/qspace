function C=minus(A,B)
% overloading - operator
% Wb,Nov17,11

  if isa(A,'SymOp') && isa(B,'SymOp')
     as=regexprep(A.istr,'^[^=]+=[ ]+','');
     bs=regexprep(B.istr,'^[^=]+=[ ]+','');

     if isempty(A), C=-B; elseif isempty(B), C=A;
     else
        C=A; C.istr=[ as ' - ' bs ];
        C.op=C.op-B.op;
        C.hc=C.hc-B.hc;
     end
  elseif isnumeric(B) && numel(B)==1 && isa(A,'SymOp')
     if B
        C=A; q=C.op;
        C.op=q-diag(repmat(B,1,size(q,1)));
        e=abs(trace(C.op))/norm(C.op); if e<1E-12
             C.istr=['rmtrace( ' C.istr ' )'];
        else C.istr=[C.istr sprintf(' %+.4g',-B)];
        end
     end
  elseif isnumeric(B) && isequal(size(A.op),size(B))
     as=regexprep(A.istr,'^[^=]+=[ ]+','');
     C=A; C.istr=[ as ' - <M>' ];
     C.op=A.op-B;
  else error('Wb:ERR','\n   ERR invalid usage'); end

end

