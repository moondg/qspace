function C=mrdivide(A,B)
% function C=mrdivide(A,B)
%
%    This overloads the right-devide operator.
%    It is restricted to numbes B only, and hence permits
%    the syntac A/x for QSpace A and x a numeric value.
%
% Wb,Sep28,12

  if nargin~=2 || ~isnumeric(B)
     helpthis, if nargin || nargout
     wbdie('invalid usage'), end, return
  end

  C=A; fac=1/B;

  for k=1:numel(C), n=length(C(k).data);
      data=C(k).data; for i=1:n, data{i}=fac*data{i}; end
      C(k).data=data;
  end

end

