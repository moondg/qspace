function tst_SdotS(ops,n)
% function tst_SdotS(ops,n)
%
%    Double check AKLT-type spin operators,
%    assuming input [S, Sa=S^n, Sb=S^n (flipped)].
%
% See also Hamilton1D/private/get_aklt_ops.m
% Wb,Sep29,15

  if isfield(ops,'op')
       S=ops(1).op; Sa=ops(2).op; Sb=ops(3).op;
  else S=ops(1); Sa=ops(2); Sb=ops(3);
  end

  if isempty(S.info.otype), S.info.otype='operator'; end

  A=permute(QSpace(getIdentityQS(S,S)),[1 3 2]);
  X=A;
  for i=1:n
     X=QSpace(contractQS(contractQS(S,'1*',X,1),'42',S,'23'));
  end
  H1=QSpace(contractQS(A,'13*',X,'13'));

  X=QSpace(contractQS(contractQS(Sa,'1*',A,1),'42',Sb,'23'));
  H2=QSpace(contractQS(A,'13*',X,'13'));

  X=QSpace(contractQS(contractQS(Sb,'1*',A,1),'42',Sa,'23'));
  H3=QSpace(contractQS(A,'13*',X,'13'));

  e=[ norm(H1-H2), norm(H1-H3) ] / dim(A,2,'-f'); % / dim(H1,1,'-f');

  if any(e>1E-12)
       wbdie('inconsistent (S.S)^%g !? (e=%.3g)',n,norm(e));
  else wblog('ok.','(S.S)^%g consistent @ e=%.3g',n,norm(e)); 
  end

return

  if norm(Sa-Sb)>1E-12
     X=QSpace(contractQS(contractQS(Sa,'1*',A,1),'42',Sa,'23'));
     H4=QSpace(contractQS(A,'13*',X,'13'));

     X=QSpace(contractQS(contractQS(Sb,'1*',A,1),'42',Sb,'23'));
     H5=QSpace(contractQS(A,'13*',X,'13'));

     wblog('TST',''); e=[ norm(H1-H4), norm(H1-H5) ]
  else
  end

end

