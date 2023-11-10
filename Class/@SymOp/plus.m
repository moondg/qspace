function C=plus(A,B)
% overloading + operator
% Wb,Nov17,11

  if isa(A,'SymOp') && isa(B,'SymOp')
     if isempty(A), C=B; elseif isempty(B), C=A;
     else C=A;
        for i=1:numel(A)
           if     isempty(A(i).op) && isempty(A(i).istr), C(i)=B(i);
           elseif isempty(B(i).op) && isempty(B(i).istr), C(i)=A(i);
           else
              as=regexprep(A(i).istr,'^[^=]+=[ ]+','');
              bs=regexprep(B(i).istr,'^[^=]+=[ ]+','');

              na=nnz(A(i).op); nb=nnz(B(i).op);
              if na && nb
                 if isempty(regexp(bs,'^[+-]'))
                      C(i).istr=[ as ' + ' bs ];
                 else C(i).istr=[ as ' ' bs ];
                 end
              elseif na, C(i).istr=as;
              elseif nb, C(i).istr=bs; end

              C(i).op=C(i).op+B(i).op;
              C(i).hc=C(i).hc+B(i).hc;
           end
        end
     end
  elseif isnumeric(B) && numel(B)==1 && isa(A,'SymOp')

     C=A; q=C.op; C.istr=sprintf('%s %+g',C.istr,B);
     C.op=q+diag(repmat(B,1,size(q,1)));

  else wbdie('invalid usage'); end

end

