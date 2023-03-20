function C=mrdivide(A,B)
% overloading * operator
% Wb,Sep28,12

  if nargin~=2 || ~isnumeric(B)
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  C=A; fac=1/B;

  for k=1:numel(C), n=length(C(k).data);
      data=C(k).data; for i=1:n, data{i}=fac*data{i}; end
      C(k).data=data;
  end

end

