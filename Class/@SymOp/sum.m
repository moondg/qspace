function C=sum(A,dim)
% overloading sum()
% Wb,Nov18,11

  if isempty(A), C=A; return; end

  s=size(A); r=numel(s);
  if nargin<2
     if r>2, error('Wb:ERR',...
     '\n   ERR invalid usage (got rank-%g object)',r); end
     if s(1)==1, A=permute(A,[2 1]); end
     Sfin=[];
  else
     A=reshape(permute(A,[dim,1:dim-1,dim+1:r]),s(dim),[]);
     Sfin=s; if dim>2, Sfin(dim)=[]; else Sfin(dim)=1; end
  end

  [n,m]=size(A);

  for j=1:m, q=A(1,j);
     for i=2:n
        q.istr=[ q.istr ' + ' A(i,j).istr ];
        q.op=q.op+A(i,j).op;
        q.hc=q.hc+A(i,j).hc;
     end
     C(j)=q;
  end

  if ~isempty(Sfin), C=reshape(C,Sfin); end

end

