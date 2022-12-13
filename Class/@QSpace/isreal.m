function isr=isreal(A)
% function i=isreal(A)
%
%    checks whether A contains complex data
%    (this overloads matlabs isreal() function).
%
% Wb,Apr23,15

% NB! use isreal rather than iscomplex, since isreal
% also corresponds to a standard matlab function.

  isr=true(size(A)); n=numel(A);

  for k=1:n, q=A(k).data;
     for i=1:numel(q)
        if ~isreal(q{i}), isr(k)=0; break; end
     end
  end

end

