function C=mtimes(A,B)
% overloading * operator
% Wb,Nov17,11

  if isa(A,'SymOp') && isa(B,'SymOp')
     C=A; C.istr=[ get_mstr(A) get_mstr(B) ];
     C.op=C.op*B.op;
     C.hc=C.hc*B.hc;
  elseif isa(B,'SymOp') && isnumeric(A), C=B;
     s=get_facstr(A,inputname(1));
     for i=1:numel(C)
        C(i).istr=[ s, get_mstr(C(i)) ];
        C(i).op=A*C(i).op; if ~isempty(C(i).hc)
        C(i).hc=A*C(i).hc; end
     end
  elseif isa(A,'SymOp') && isnumeric(B), C=A;
     s=get_facstr(B,inputname(2),'R');
     for i=1:numel(C)
        if numel(B)==1
             C(i).istr=[ s, get_mstr(C(i))    ];
        else C(i).istr=[    get_mstr(C(i)), s ];
        end
        C(i).op=C(i).op*B; if ~isempty(C(i).hc)
        C(i).hc=C(i).hc*B; end
     end
  else wbdie('invalid usage'); end

end

function s=get_facstr(x,s,rflag)

  if nargin<3, rflag=0; end

  if 1 || isempty(s)
     if numel(x)==1
        if     x==+1, s='';
        elseif x==-1, s='-';
        else          s=sprintf('%.4g ',full(x)); end
     else
        s=sprintf('x%g',size(A)); s=s(2:end);
        if rflag, s=['*(' s ')']; else s=['(' s ')*']; end
     end
  else
     if rflag, s=[' * ',s]; else s=[s,' * ']; end
  end

end

