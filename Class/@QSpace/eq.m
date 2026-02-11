function i=eq(A,B)
% function i=eq(A,B)
%
%    Simple wrapper to QSpace/isequal.m (whether two QSpaces are equal)
%    that overloads Matlab's the == (equal) operator.
%
% Wb,Feb22,24

  i=isequal(A,B);

  if ~i && isequal(size(A),size(B)), n=numel(A); i=1;
    for k=1:n
       if ~isequal(A(k).Q,   B(k).Q   ), i=0; break; end
       if ~isequal(A(k).data,B(k).data), i=0; break; end
    end
    if ~i, return; end

    for k=1:n
       Ia=A(k).info; if isfield(Ia,'ctime')
       Ib=B(k).info; Ib.ctime=Ia.ctime;
       if ~isequal(Ia,Ib), i=0; break; end; end
    end
  end

end

