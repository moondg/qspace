function kc=initPsi(HAM)
% function kc=initPsi(HAM)
%
%    initialize A to random state with very few symmetry
%    sectors along the way (close to product state).
%
%    If HAM is requested as output argument,
%    the state A will be brought into the local form of kc=1,
%    with the matrix represenation for the Hamiltonians
%    generated (i.e. the routine initHam() is called).
%    (HAM is altered only in that kc is updated!).
%
% Wb,Apr08,14

  L=numel(HAM.mpo);
  A=QSpace(1,L);

  E=HAM.ops(1);
  if ~isIdentityQS(E)
     wbdie('invalid HAM (expecting');
  end

  A(1)=addsingleton(E,'-p');
  for k=2:L
     Q=QSpace(getIdentityQS(A(k-1),2,E));
     d3=getzdim(Q,3,'-p') + 1E-6*sqrt(sum(Q.Q{3}.^2,2));
     [d3,is]=sort(d3); if numel(is)>3, is=is(1:3); end
     if k==L, is=is(1); end

     Q=getsub(Q,is);
     for i=1:numel(Q.data)
        x=Q.data{i}(:,:,1); x=randn(size(x));
        Q.data{i}=x/norm(x);
     end

     R=contractQS(Q,[1 2],Q,[1 2]);
     for i=1:numel(R.data)
        [u,r]=eig(R.data{i}); r=diag(r); j=find(r>1E-8);
        R.data{i}=u(:,j)*diag(1./sqrt(r(j))); % *u';
     end

     A(k)=contractQS(Q,3,R,1,[1 3 2]);
  end

  A=A(L:-1:1);
  for k=1:L
     t3=sprintf('K%02g,K%02g*,s%02g',k-1,k,k);
     A(k)=permuteQS(A(k),[2 1 3]);
     A(k).info.itags=strread(t3,'%s','whitespace',',')';

     Xk=load_dmrg_data(HAM,k);
        Xk.AK=A(k);
        Xk.info.odir=1;
     save_dmrg_data(HAM,Xk,k);
  end

  if nargout, kc=1; end

end

