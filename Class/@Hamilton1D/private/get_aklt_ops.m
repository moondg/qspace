function Sab=get_aklt_ops(S,naklt,waklt,sstr)
% function Sab=get_aklt_ops(S,naklt,waklt,sstr)
%
%    get paired-up operators Sa(2*k+1).Sa(2*k) for (S.S)^(k+2)
%    having k=0:naklt-1.
%
%    In order to deal with OM in the irop space on a generic level,
%    fully contract S'.S term, then split it again to obtain minimal
%    irop space.
%
% DEPRECATED (see get_poly_ops) // Wb,Apr02,20
%
% Wb,Sep25,15 ; Wb,Feb19,16

% for old version that combines irop indices for the operators
% S^(k+2) individualy, see Archive/get_aklt_ops_160219.m
% outsourced from setup_Heisenberg.m // Wb,Sep25,15

  if nargin<2 || naklt<1 || nargout>1
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end
  if nargin<3, sstr='S'; end

% -------------------------------------------------------------- %
% NB! AKLT terms (as well as Heisenberg interactions between     %
% different sites with different spin S require the contraction  %
% Sa.Sb' with Sa!=Sb (!)  // Wb,Aug27,15                         %
% -------------------------------------------------------------- %

  SS=permuteQS(contractQS(S,'3*',S,3),'2314');
  SX=SS;

  Z1=getIdentityQS(SS,3,'-0');
  A2=getIdentityQS(SS,1,Z1);
  A2=QSpace(contractQS(A2,2,Z1,'2*',[1 3 2]));

  for ik=1:naklt
     SX=QSpace(contractQS(SS,'34',SX,'12'));

     X2=QSpace(contractQS(contractQS(A2,'12*',SX,'13'),'23',A2,'12*'));

     [u,s,v]=svdQS(X2,2,'stol',1E-12);

     e=getIdentityQS(s,1); if ~isequal(e.Q{1},s.Q{1})
        wbdie('got different Q sorting !?'); end

     for i=1:numel(s.data)
        si=s.data{i}; if any(si)<=0
           wbdie('got min(svd)=%.3g !?',min(si)); end
        s.data{i}=diag(sqrt(si));
        e.data{i}=s.data{i};
     end

     Sa=skipzeros(QSpace(contractQS(A2,3,contractQS(u,2,e,1),1)));
     Sb=skipzeros(QSpace(contractQS(A2,3,contractQS(s,2,v,1),2)));

     Sb=QSpace(permuteQS(Sb,'213*'));

     if 0
        Q=QSpace(contractQS(Sa,'3*',Sb,'3',[2 3 1 4])) - SX;
        e=normQS(Q)/normQS(SX); wblog('-->','e=%.5g',e); 
     end

     Sab(1,ik)=init_ops(Sa,sprintf('spin %s^%g [a]',sstr,ik+1),'~hconj');
     Sab(2,ik)=init_ops(Sb,sprintf('spin %s^%g [b]',sstr,ik+1),'~hconj');

     tst_SdotS([S,Sab(1,ik).op,Sab(2,ik).op],ik+1);
  end

  Sab=Sab(:);

end

% -------------------------------------------------------------------- %

