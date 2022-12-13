function A=compact(A)
% merge quantum numbers of tensor operator
% Wb,Jul07,09

  if isempty(A.Q), return, end

  r=length(A.Q);
  if mod(r,2), error('Wb:ERR','invalid usage (rank=%g)',r); end

wblog('WRN','double check this');

  B=A+A';
  for i=1:r, qq{i}=expandQ(B,i); end

  i1=1:(r/2); i2=(r/2)+1:r;

  E1=QSpace(qq{i1},'-Rlast','identity');
  E2=QSpace(qq{i2},'-Rlast','identity');

  A=QSpace(contractQS(E1,i1,contractQS(A,i2,E2,i1),i1));

end

