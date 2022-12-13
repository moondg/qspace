function S=sum(A,dim)
% overloading the sum operator

  if isempty(A), S=A; return; end

  if nargin<2 && isvector(A)
     S=struct(A(1));
     for i=2:length(A), S=plusQS(S,A(i)); end
  else
     if nargin<2, dim=1;
     elseif dim~=1 && dim~=2, dim, error('Wb:ERR','invalid usage'); end

     [n,m]=size(A);

     if dim==1
        S=struct(A(1,:));
        for j=1:m, for i=2:n, S(j)=plusQS(S(j),A(i,j)); end, end
     else
        S=struct(A(:,1));
        for i=1:n, for j=2:m, S(i)=plusQS(S(i),A(i,j)); end, end
     end
  end

  S=skipzeros(class(S,'QSpace'));

end

