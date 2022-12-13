function C=plus(A,B)
% overloading the + operator

  if isnumeric(B) && isscalar(B)
   % adding single scalar = acts like diagonal operator
     if isempty(A), C=A; return; end
     if length(A.Q)~=2, error('Wb:ERR',...
       '\n   ERR invalid usage (got rank-%d QSpace)',numel(A.Q)); end

     [isd,s]=isdiag(A);
     if ~isempty(s), error('Wb:ERR',estr); end

     C=A; n=length(C.data);

     if isd>1
        for i=1:n, C.data{i} = C.data{i} + B; end
     else
        for i=1:n
           if isequal(C.Q{1}(i,:),C.Q{2}(i,:))
              C.data{i} = C.data{i} + full(B*speye(size(C.data{i})));
           end
        end
     end
  else
     if ~isequal(size(A),size(B)), error('Wb:ERR',...
        '\n   ERR invalid usage (size mismatch between A and B)'); end
     C=A; n=numel(A);
     for i=1:n, C(i)=QSpace(plusQS(A(i),B(i))); end
  end

end

