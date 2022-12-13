function A=skiprecs(A,ix)
% function A=skiprecs(A,ix)
% Wb,Aug01,12

% adapted from getrecs()

  for p=1:length(A.Q), A.Q{p}(ix,:)=[]; end

% reduce cgr first (since gotCGS checks length(data))
  if gotCGS(A), A.info.cgr(ix,:)=[]; end

  A.data(ix)=[];

end

