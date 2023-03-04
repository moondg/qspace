% Wb,Mar03,23

  [FF,Z,SS,IS]=getLocalSpace('FermionS','Acharge,SU2spin,SU3channel','NC',3);
  E = IS.E

  Lsep=repmat('─',1,90);
  fprintf(1,'%s\n\n',Lsep); % ------------------------------------------------

  A=getIdentity(E,E);
  A=setitags(A,{'K1','s2','K2'});

  display(A,'-c')
  fprintf(1,'%s\n\n',Lsep); % ------------------------------------------------

  H=contract(A,'!3*',{SS,'-op:K1','*',{A,SS,'-op:^s'}})
  fprintf(1,'%s\n\n',Lsep); % ------------------------------------------------

  i=find(H.Q{1}(:,2)==6); getsub(H,i)
  fprintf(1,'%s\n\n',Lsep); % ------------------------------------------------

  if H.data{i} == 9/4
     wblog('ok!','got correct S.S interaction @ %g :)\N',H.data{i});
  end

