% Setup for chiral two-impurity problem coupled to two `bath channels' %
% which also permits `cross couplings' due to off-diagonaly            %
% hybridization, so they effectively describe the same bath.           %
% Therefore no channel symmetry is used, using global charge and spin  %
% symmetries only // Wb,Feb01,21                                       %
% -------------------------------------------------------------------- %
% adapted from setupSIAM_SU2x2.m

  global param

  setdef('U',0.08, 'Gamma',0.01, 'Lambda',4, 'N',30);
  setdef('epsd',-U/2, 'B',0);

  istr=''; sym='';
  param=add2struct('-',istr,sym,U,epsd,Gamma,B,Lambda,N);

  if exist('TK0')==1, TK=TK0; else
     TK=TKondo;
  end

  param.istr='two-channel chiral Anderson model';

% -------------------------------------------------------------------- %
% operator setup
% -------------------------------------------------------------------- %
% this is where the symmetry is set once and for all:
% the largest symmetry is given for SU(2) spin symmetry (B=0)
% and SU(2) particle-hole (p/h) symmetry (epsd=-U/2) for each impurity
% while *not* assuming channel symmetry (this is therefore also
% the fastest 2-channel calculation on practical grounds, minutes
% vs. tens of minutes)
%
%  * Because no channel symmetr i asssumed
%    Gamma and U can differ for each impurity
%
%  * When swichting to Kondo setup with spin-spin interactions
%    for each impurity, their respective Kondo coupling can also
%    be chosen differently. When simulating the Kondo Hamiltonian
%    itself, the coupling terms to the 0-th bath sites need to be
%    included with H0 itself (this is options for the Anderson-
%    type impurities below, hence not yet done; nevertheless,
%    the Kondo setting with isotropic S.S interactions can be
%    simply mimicked in a simple manner by taking U, Gamma >> D
%    (with D=1 the half-bandwidth) with J_Kondo^effective ~ Gamma/U
%    in the Anderson model itself.
%
% -------------------------------------------------------------------- %

  setdef('SYM',''); clear ZFLAG*

  os={'FermionS','','NC',2,'-v'};

  if epsd~=-U/2 && B~=0 || isequal(SYM,'A,A')
     os{2}='Aspin,Acharge';
     [FF,Z,SS,IS]=getLocalSpace(os{:});

     NF=FF; for i=1:numel(FF) 
        NF(i)=contractQS(FF(i),'13*',FF(i),'13');
     end

     HU = NF(1)*NF(3) + NF(2)*NF(4);

     H0 = epsd*sum(NF) + U*HU - B*SS(3);

  elseif epsd~=-U/2 && B==0 || isequal(SYM,'SU2,A')
     os{2}='SU2spin,Acharge';
     [FF,Z,SS,IS]=getLocalSpace(os{:});

     NF=FF; for i=1:numel(FF) 
        NF(i)=contractQS(FF(i),'13*',FF(i),'13');
     end

     HU=0.5*(NF(1)*NF(1)-NF(1)) ...
      + 0.5*(NF(2)*NF(2)-NF(2));

     H0 = epsd*sum(NF) + U*HU;

  elseif epsd==-U/2 && B==0 || isequal(SYM,'SU2,SU2')
     os{2}='SU2spin,SU2charge';
     [FF,Z,SS,IS]=getLocalSpace(os{:}); ZFLAG=2;

     HU=(1/2)*[(4/3)*sum(IS.Q2)];

     H0 = U*HU;

  else wberr('invalid usage');
  end

  E0=IS.E;

  if isempty(H0), H0=0*E0; end

  q=getvac(Z);
  A0=QSpace(permuteQS(getIdentityQS(q,1,Z,1),[1 3 2]));

  SYM=getsym(H0); param.sym=IS.sym;

  co={'-x'};
  if exist('z','var'), co=[co {'z',z}]; param.z=z; end
  if exist('ALambda','var') && ~ALambda
     co(end+1)={'-w'};
     if isfield(param,'ALambda'), param=rmfield(param,'ALambda'); end
  end

  [ff,~,Iwc]=getNRGcoupling(Gamma,Lambda,N,co{:});

  ff=ff(:)*reshape(eye(numel(FF)),1,[]);
  % ff as a single vector, only uses diagonal couplings,  i.e.
  % equivalent to ff(k,i,i) above; only non-zero couplings ff
  % are explicitly contracted within NRGWilsonQS.
  %
  % extending to a matrix above permits one to specify
  % the precise hopping terms ff(k,i,j) coupling FF(i)
  % at Wilson shell k-1 to FF(j) in the new shell k;
  % (fdmNRG_QS requires ff n matrix format, hence the
  % reshape into ff(:,:) which fuses (i,j) with i the
  % fast index; above by using eye(), this is equivaletn
  % matter of fact, to just haveing left ff as the
  % initial vector.

% -------------------------------------------------------------------- %
  FC=FF;
  Z0=Z;

  op1=[];
  op2=[ FC(1), SS(end) ];

  zflags=[1 0];

% -------------------------------------------------------------------- %

  inl(1); disp(param);

