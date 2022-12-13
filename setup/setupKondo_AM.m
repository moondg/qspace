% = README =========================================================== %
%
% * NC=1: J>2*Gamma makes fix point go to infinity VERY quickly
%
% * sign of J (checked for NC=1)
%   J<0 takes spinful state as groundstate : very sensitive to B!=0
%   J>0 takes spinless state (singlet/double occupancy) as groundstate
%       making the system rather insensitive to B!=0 (checked)
%
% * NC=2: J=-0.02, Gamma=0.01 (same as for NC=1) has much smaller TK!
%       checked for Nkeep=4096 (takes 42GB of harddrive -> see frankie)
%       not converged for Nkeep=512; Nkeep=1024 appears close to
%       sufficient (i.e. Eflow nearly the same as for Nkeep=4096)
%
% * estimate of how many states to keep
%   -> NC=1: Nkeep=256 appears close to sufficient
%   -> NC=2: 256*4=1024 required at minimum (seems to agree with above)
%   -> NC=3: 256*4*4=4096 required at minimum: seems doable
%
% Wb,Jul26,07
% ==================================================================== %

% -------------------------------------------------------------------- %
% setup operator set for multi-channel Kondo problem
% in Anderson Hamiltonina language
% [ Discussion with A Rosch, Jun08,07 ]
% -------------------------------------------------------------------- %

  global param f4

  if ~exist('J','var') || ~exist('Gamma','var') || ~exist('Lambda','var')
     wblog('TST Setting test parameters.');
     U=0; Gamma=0.01; epsd=-U/2; B=0E-6; N=50; J=-0.02; Lambda=2;

     NC=2; Nkeep=512;
  end

  param=add2struct('-',U,epsd,Gamma,J,B,Lambda,N);

  co={};
  if exist('z',   'var'), co(end+1:end+2)={'z',z   }; param.z=z;    end
  if exist('pg_r','var'), co(end+1:end+2)={'r',pg_r}; param.r=pg_r; end
  if exist('wflag','var') && wflag, co={co{:},'-w'}; end

  ff=getNRGcoupling(Gamma,Lambda,N,co{:});

  if exist('TK0')==1, TK=TK0; else
  TK=TKondo; end

% ==================================================================== %
  initdef('NC',1);

  wblog('<i>', 'Kondo parameters (NC=%g, TK=%.4g)\N',NC,TK);
  disp(param);

  sx=spinmat('-pauli');
  sx{2}=-1i*sx{2};

  QS=[  0   0
        1  +1
        1  -1
        2   0 ];

  Qloc=repmat({zeros(4,NC+1)},1,NC);
  for i=1:NC, Qloc{i}(:,[i end])=QS; end

  f4.c=QSpace(2,NC); f4.cn=QSpace(2,NC);
  f4.z=QSpace(1,NC); f4.sx=QSpace(3,NC);

  for l=1:NC, Ql=Qloc{l};
     f4.z(l)=QSpace(Ql, diag((-1).^[0 1 1 2]), 'operator');

     f4.c(:,l)=[
         QSpace(...
            Ql([1 2],:),  1, ...
            Ql([3 4],:),  1  ...
         )
         QSpace(...
            Ql([1 3],:),  1, ...
            Ql([2 4],:), -1  ...
         )
     ];

     f4.cn(:,l)=[
         QSpace(Ql([3 4],:),  1)
         QSpace(Ql([2 4],:), -1)
     ];

     f4.sx(:,l)=[
        QSpace(Ql(2:3,:),sx{1},'operator')
        QSpace(Ql(2:3,:),sx{2},'operator')
        QSpace(Ql(2:3,:),sx{3},'operator')
     ];
  end

% -------------------------------------------------------------------- %
  [A0,Aloc]=QSpace(zeros(1,NC+1), Qloc{:}, 'identity');

  A0=QSpace(contractQS(A0,3:2+NC,Aloc,1:NC));

  E4=QSpace(contractQS(Aloc,1:NC,Aloc,1:NC));

  FC=QSpace(2,NC); SX=QSpace(3,NC); Q=Aloc;
  for l=1:NC,
     iperm=cmat_iperm(NC+1,l);
     ix=iperm(1:end-1);

     for i=1:2
         FC(i,l)=...
            contractQS(Aloc, 1:NC,...
            contractQS(Q,l,f4.c(i,l),2), ix,'conjA'...
         );
     end
     Q=contractQS(Q,l,f4.z(l),2,iperm);

     for p=1:3
         SX(p,l)=...
            contractQS(Aloc, 1:NC,...
            contractQS(Aloc,l,f4.sx(p,l),2), ix,'conjA'...
         );
     end
  end

  Z=QSpace(contractQS(Aloc, 1:NC, Q, 1:NC,'conjA'));

  H0=0*E4; if exist('B','var') && B~=0
  for l=1:NC, H0=H0+(B/2)*SX(3,l); end, end

  SS=QSpace;
  for p=1:3, Ql=QSpace; for l=1:NC, Ql=Ql+SX(p,l); end
     if p~=2
          SS=SS+Ql*Ql;
     else SS=SS-Ql*Ql; end
  end
  H0=H0+J*SS;

% -------------------------------------------------------------------- %
  Z0=Z;

  i=1;
  op1=[FC(i); comm(FC(i),SS)];
  op2= FC([i;i]);

  SZ=sum(SX(3,:));

  if exist('chiflag','var') && chiflag
     op1(3)=SZ; zflags=[1 1 0];
     op2(3)=SZ; cflags=[0 0 1];
  end

  FC_=FC;

% -------------------------------------------------------------------- %
  clear i k l p ic ix f m2 Q data co iperm Ql
  e=0;
  for l=1:NC, for k=l:NC
  for i=1:2,  for j=i:2
     a=[ acomm(FC(i,l),FC(j,k)), acomm(FC(i,l),FC(j,k)') ];
     if i==j && l==k
        if ~isIdentityQS(a(2))
        wblog('ERR','invalid CR for l=%g',l); e=e+1; end
     else
        if normQS(a(2))~=0
        wblog('ERR','invalid CR for l=%g',l); e=e+1; end
     end
        if normQS(a(1))~=0
        wblog('ERR','invalid CR for l=%g',l); e=e+1; end
  end, end
  end, end

  for l=1:NC, for k=l:NC
  for i=1:3,  for j=i+1:3
     a=comm(SX(i,l),SX(j,k));
     if l~=k
        if normQS(a)~=0
        wblog('ERR','invalid CR [%s]',vec2str([i l j k])); e=e+1; end
        continue
     end

     i3=1:3; i3([i j])=[];

     if normQS(a-2*SX(i3,l))~=0
     wblog('ERR','invalid CR [%s]',vec2str([i l j k])); e=e+1; end

  end, end
  end, end

  if e==0, wblog('SUC',':)'); end

  clear l i j k e a

