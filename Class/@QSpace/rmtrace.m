function A=rmtrace(A)
% function A=rmtrace(A)
%
%    remove trace from input tensor A
%
% Wb,Apr24,13

  for i=1:numel(A)
     if numel(A(i).Q)~=2, error('Wb:ERR',...
       '\n   ERR invalid usage (rank-2 tensors only)'); end
     A(i)=A(i)-trace(A(i))/dim(A(i),1,'-f');
  end

end

