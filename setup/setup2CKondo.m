
% -------------------------------------------------------------------- %
% operator setup for two-channel Kondo problem
% NB! see also setupKondo with NC=2 (!)
% NB! see also setupKondo_SU2x2.m for 1-channel setup
% Wb,Mar01,08
% Wb,Apr24,12 : adapted to non-abelian symmetries [SU(2)_spin * Sp(4)]
% -------------------------------------------------------------------- %

% NC=2; % two channels (modified from setupKondo_AM)

% NB! J=0.12 is TOO SMALL (TK<1E-8)
  setdef(...
     'J', 0.25, 'B',0E-6, 'NC',2, ...
     'N',36, 'Lambda',4,'Nkeep',1024,'Etrunc',8 ...
  );

  if numel(J)>1
       dJ=diff(J); if ~dJ, J=J(1); end
  else dJ=0; end

  global param
  param=add2struct('-',J,B,Lambda,N,Nkeep,Etrunc);
  TK=TKondo;

  if exist('z','var'), co={'z',z}; param.z=z; else co={}; end
  ff=getNRGcoupling(1,Lambda,N,'-x',co{:});
  ff=ff(2:end);

  wblog('<i>', 'Kondo parameters (NC=%g, TK=%.4g)\N',NC,TK);

  setdef('SYM','');
  if B || ~isempty(regexp(SYM,'^A,'));
       sym={'Aspin',   'SpNchannel'};
  elseif dJ, sym={'SU2spin', 'SU2charge(:)'};
  else sym={'SU2spin', 'SpNchannel'};
  end

  if ~isempty(SYM)
     if     ~isempty(regexp(SYM,',Sp'))       sym{2}='SpNchannel';
     elseif ~isempty(regexp(SYM,',A,SU'))     sym{2}='Acharge,SUNchannel';
     elseif ~isempty(regexp(SYM,',A,A'))      sym{2}='Acharge(:)';
     elseif ~isempty(regexp(SYM,',SU2,SU2'))  sym{2}='SU2charge(:)';
     end
  end

  [~,~, s0, ~]=getLocalSpace('FermionS', sym{1},'NC',1);
  [FC,Z,SS,IS]=getLocalSpace('FermionS',[sym{1} ',' sym{2}],'NC',NC,'-v');

  SYM=getsym(IS.E); param.sym=IS.sym;

  structdisp(param);

  if ~isempty(regexp(sym{2},'^(Sp|SU)'))
       ZFLAG=3;
  else ZFLAG=1; end

  s0_=s0;

  s=getsym(FC(1),'-c');
  for i=2:numel(s)
     s0=addSymmetry(s0,s{i});
  end

  A0=getIdentity(s0(end),Z,[1 3 2]);

  r=numel(s0);
  if r~=1 && r~=3, wbdie('invalid spin setting'); end

  if dJ, SS_=SS; SS=IS.S3; end

  [n,m]=size(SS); q=[ numel(J), n ];
  if diff(q), wbdie('invalid J (len=%d / %d)',q); end
  if m~=r, wbdie('invalid spin setting'); end

  H0=QSpace;
  for p=1:n
  for i=1:m, Q=contract(J(p)*s0(i),'1*',A0,1);
     if numel(s0(i).Q)==3
          Q=contract(Q,'42',SS(p,i),'23');
     else Q=contract(Q,3,SS(p,i),2);
     end
     Q=contract(A0,'13*',Q,'13');
     H0=H0+Q;
  end
  end
  HJ=H0;

  if B
     H0=H0+contract(A0,'13*',contractQS(-B*s0(end),2,A0,1),'13');
  end

  FX=FC; OX=FC;
  for i=1:numel(FC)
     FX(i)=contract(A0,'13*',contractQS(A0,3,FC(i),2),'13');
     OX(i)=QSpace(contractQS(FX(i),2,HJ,1,[1 3 2])) - contractQS(HJ,2,FX(i),1);
  end

  op1=[];
  op2=OX; op2(end+1)=FX(1);

  op2(end+1)=contractQS(A0,'13*',contractQS(A0,1,s0(end),2),'32');

  for i=1:numel(op2)
     if numel(op2(i).Q)>2
     op2(i).info.otype='operator'; end
  end
  zflags=ones(1,numel(op2)); zflags(end)=0;

  Z0=contract(A0,'13*',contractQS(A0,3,Z,2),'13');

  AJ=A0; H0_=H0;
  A0=getIdentity(Z0,Z,[1 3 2]);

  f1=ff(1); ff=ff(2:end);

  H0=contract(A0,'13*',contractQS(H0,2,A0,1),'13');

  for i=1:numel(FC)
     if ZFLAG>1
        H0=H0+contract(A0,'13*',...
           contractQS(contractQS(FX(i),'1*',A0,1),'42',...
           f1*FC(i),'23'), ...
        '13');
     else
        Q=contract(A0,'13*',...
           contractQS(contractQS(FX(i),'1*',A0,1),'42',...
           contractQS(f1*Z,2,FC(i),1),[2 3]), ...
        '13');
        H0=H0+Q+Q';
     end
  end

