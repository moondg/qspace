function dd=catdiag2(varargin)
% function dd=catdiag2(A,B,...)
%
%    concatenate all diag(A(i)) as columns of a matrix
%    i.e. along dim=2.
%
% Wb,Dec05,11

  for k=1:nargin
     A=varargin{k}; n=numel(A); d=cell(n);
     for i=1:n, d{i}=diag(A(i).op); end
     dd{k}=cat(2,d{:});
  end

  if nargin>1, dd=cat(2,dd{:}); else dd=dd{1}; end

end

