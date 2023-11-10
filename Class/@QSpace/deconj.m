function [A,t]=deconj(A,idx)
% Function [A,t]=deconj(A [,idx])
%
%    Unset "bra" (complex conjugate) label for specified dimensions
%    (all if idx is not specified).
%
% Wb,Apr17,13

  if nargin<2, idx=1:numel(A(1).Q); end
  n=numel(idx);

  c='*';

  for k=1:numel(A)
     [i,t]=gotITags(A(k)); m=0;
     if n>i || max(idx)>i
        wbdie('input index set out of bounds or not unique');
     end
     for j=1:n, i=idx(j);
        if t{i}(end)==c, t{i}=t{i}(1:end-1); m=m+1; end
     end
     if m, A(k).info.itags=t; end
  end

end

