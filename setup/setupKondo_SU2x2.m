% setup operator set for standard Kondo problem (single channel)
% using particle-hole and spin SU(2) symmetry (if B=0)
% Wb,Sep19,11
% -------------------------------------------------------------------- %

% adapted from setupSIAM.m

  if ~isvar('J') || isempty(J)
     if isvar('U') && isvar('Gamma')
          J=8*Gamma/(pi*U);
     else J=0.12; end
  end
  setdef('B',0,'Lambda',2,'N',50);

  istr=''; sym=''; global param
  param=add2struct('-',istr,sym,'fJ=J',B,Lambda,N);

  co={'-x'};
  if exist('z',   'var'), co(end+1:end+2)={'z',z   }; param.z=z;    end
  if exist('pg_r','var'), co(end+1:end+2)={'r',pg_r}; param.r=pg_r; end
  setopts(co,'Nmax?');

  if exist('ALambda','var')
     if ischar(ALambda) || isempty(ALambda) || ALambda
        co(end+1)={'-AL'}; param.ALambda=1;
     else
        co(end+1)={'-w'};
        if isfield(param,'ALambda') param=rmfield(param,'ALambda'); end
     end
  end

  [ff,q,Ic]=getNRGcoupling(1,Lambda,N,co{:});
  ff(1)=J;

  if exist('TK0')==1, TK=TK0; else
  TK=TKondo; end

% -------------------------------------------------------------------- %
  setdef('NC',1);
  param.NC=NC; o={'NC',NC,'-v'};

  param.istr=sprintf('%s-channel Kondo(J) model',num2text(NC));

  Bflag=(isset('B') || (isvar('Jz') && Jz~=J));

  if NC>1
     if Bflag, sp='Aspin'; else sp='SU2spin'; end
     if 1
          sym=[           sp ',SpNchannel'];
     else sym=['Acharge,' sp ',SUNchannel'];
     end; clear sp

     [FF,Z,SS,IS]=getLocalSpace('FermionS',sym,o{:});

  elseif isset('SYM') && isequal(SYM,'A,A')
     [FF,Z,SS,IS]=getLocalSpace('FermionS','Acharge,Aspin',  o{:});
  elseif isset('SYM') && isequal(SYM,'SU2,A') || Bflag
     [FF,Z,SS,IS]=getLocalSpace('FermionS','SU2charge,Aspin',  o{:});
  elseif isset('SYM') && isequal(SYM,'A,SU2')
     wbdie('invalid SYM=%s for Kondo model !??',SYM);
  else
     [FF,Z,SS,IS]=getLocalSpace('FermionS','SU2charge,SU2spin',o{:});
  end

  SYM=getsym(FF); param.sym=IS.sym;

  E1=QSpace(getIdentityQS(SS(end)));

  A0=QSpace(permuteQS(getIdentityQS(E1,1,Z,1),[1 3 2]));
  A1=QSpace(permuteQS(getIdentityQS(A0,2,Z,1),[1 3 2]));

  j3=J; n=numel(J);
  if n==1
     if isvar('Jz')
        j3=[J J Jz]; param.fJ=j3;
     end 
  elseif n~=3, J, wbdie('invalid usage');
  end

  j3=2*j3;

  n=numel(SS);
  if n==1
     if numel(j3)>1 && ~norm(diff(j3(:))), j3=j3(1); end
  elseif n~=3, wbdie('(double check; added Jz'); end

  if n>1, H0=QSpace;
     if numel(j3)==1, j3=repmat(j3,1,n); end
     for p=1:n, q=SS(p);
        if rank(q)>2
           q=contract(...
              contract(q,2,A0,1),'142*',...
              contract(A0,3,q,2),'134'  ...
           );
        else
           q=contract(...
              contract(q,2,A0,1),'13*',...
              contract(A0,3,q,2),'13'  ...
           );
        end
        H0=H0+j3(p)*q;
     end
  else
     if numel(j3)~=1, j3, wbdie('invalid usage'); end
     H0=j3*QSpace(contractQS(...
        contractQS(SS,2,A0,1),'142*',...
        contractQS(A0,3,SS,2),'134'  ...
     ));

     if 0 && isset('Jph'), if Jph~=1, j3=2*Jph; end
        H0=H0+j3*QSpace(contractQS(...
           contractQS(IS.C3,2,A0,1),'142*',...
           contractQS(A0,3,IS.C3,2),'134'  ...
        ));
        param.Jph=j3/2;
     end
  end

  ff=ff(2:end);

% -------------------------------------------------------------------- %
  if isset('B')
     if numel(SS)<3 || dim(SS(end),'-op')~=1
        wbdie('invalid symmetry/spin setting for B!=0');
     end
     q=SS(end); if rank(q)>2
        if dim(q,'-op')~=1
           wbdie('unexpected rank-%g spin operator',rank(q)); end
        q=fixScalarOp('-f',q);
     end
     q=QSpace(contractQS(contractQS(q,2,A0,1),'13*',A0,'13'));
     H0=H0-B*q;
  end

% -------------------------------------------------------------------- %

  op1=[]; op2=FF;
  if isset('B'), m=numel(FF); else m=1; end
  for i=1:m
     q=QSpace(contractQS(A0,'13*',contractQS(A0,3,FF(i),2),'13'));
     op2(i)=QSpace(contractQS(H0,2,q,1)) - contractQS(H0,1,q,2,[2 1 3]);
     op2(i).info.otype=FF(i).info.otype;
  end
  i=m+1;

  FC=op2;

  q=SS(end); if rank(q)>2, p={[1 3 4 2]}; else p={}; end
  op2(i)=QSpace(contractQS(A0,'13*',contractQS(q,2,A0,1,p{:}),'13'));
  op2(i).info.otype=q.info.otype;

  zflags=ones(1,m); zflags(m+1)=0;

  Z0=QSpace(contractQS(A0,'13*',contractQS(A0,3,Z,2),'13'));

% -------------------------------------------------------------------- %

  if SYM(1)=='S', ZFLAG=2; else clear ZFLAG; end
  if isequal(SYM,'A,A')
     fixAbelian(FF,FC,SS,op1,op2,'-xop3');
  end

  wblog('<i>', 'Kondo parameters (TK=%.4e)\N', TK);
  disp(param);

% -------------------------------------------------------------------- %

