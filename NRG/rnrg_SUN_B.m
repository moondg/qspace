
% same as rnrg_SUN but including magnetic field B
% => B break spin SU(2) which consequently becomes plain abelian.

% mat='Wb100918_NRG_SU3_B12_D3072_L40_J0230';

% clear ; savepara=1 ; rnrg_SUN_B
  setdef('savepara',0);

  inl(1);
  if ~exist('param','var')
     NC=3;

   % Hund's coupling: 0.008 * 2.5/(3+0.5) * 4 = 0.0229
   % see rkondo.m for NC=3;
   % NB! factor 4 since Sop is defined without factor of 0.5 there!

   % JH=0.0023;
   % JH=0.0012; % too small (looses significance compared to Gamma!!)

     JH=0.0230;
   % JH=0.0300; % TK = ~3E-8 (Lambda=4 : converged at N>30; TK ~ k=25)
   % JH=0.0500; % TK = ~3E-12 (Lambda=4 : converged at N>45; TK ~ k=38)

     Gamma=0.01; Lambda=4; L=50;
   % B=1E-12;
   % B=1E-4;

   % expected TKondo (see tex notes)
     global g_TKondo
   % g_TKondo=sqrt(JH*Gamma/2)*exp(-JH/(2*Gamma) * (1-1/(2*NC+1)));
     g_TKondo=sqrt(JH*Gamma/2)*exp(-(NC+1)*JH/Gamma); % 4=2S+1 = 2(S+1/2)
     
   % ============================= SU123 ========================== runtime  IO-time %
   % Nkeep= 256; % NB! keeps up to  5070 states including z-labels    8 min
   % Nkeep= 512; % NB! keeps up to 10050 states including z-labels   22 min  >27 min!?
     % Nkeep=512 => 3.1G HD space @ 10Mb/sec => expecting 5 min!
     % (compare previous Abelian calculation: Nk=4096: runtime 20 hrs
     % eg. see $SCR/Data/goedel rkondo_070919_sge77429.log

   % Lambda=4; setdef('Nkeep',3072);
   % Nkeep=4096;
   % Nkeep=3072;
   % Nkeep=2872;
   % Nkeep=2048;
   % Nkeep=1024;
   % Nkeep=64;   % for testing purposes

   % Lambda=3;   L=30; Nkeep=2048; % Nkeep=2560;
   % Lambda=2.5; L=30; Nkeep=2560; % Nkeep=3072;

   % test parameters (Markus Hanl), Wb,Mar08,11
     B=0; Nkeep=1000; T=3E-4;

   % ============================= SU123 ========== for comparison: (SU2^4) %
   % Nkeep=1024; % NB! keeps up to 23500 states including z-labels
   % Nkeep=1500; % NB! keeps up to 31960 states including z-labels  (23500)
   % Nkeep=2048; % NB! keeps up to 53600 states including z-labels  (32800)
   % Nkeep=2600; % NB! keeps up to 63500 states including z-labels  (40400)
   % Nkeep=3500; % NB! keeps up to 63500 states including z-labels
   % ====================================================================== %

     global param
     param=add2struct('-',Gamma,JH,'B?',Lambda,L,NC,Nkeep);

     ff=getNRGcoupling(Gamma,Lambda,L);
     param.ff=ff;
  end

% setting up SU(2)^4 for symmetric Anderson/Kondo hybrid
% Wb,Jun10,10

  f2=[0 1; 0 0]; % [1,0] = empty state
  z2=f2*f2'-f2'*f2; n2=f2'*f2; e2=eye(2);

% F-operators
  FF=cell(NC,2); n3=zeros(4^NC);
  for i=1:NC
  for j=1:2
     q=repmat({e2},1,2*NC); i1=2*(i-1)+j; q{i1}=f2;  q(i1+1:end)={z2};
     FF{i,j}=mkron(q{:});

     q=repmat({e2},1,2*NC); i1=2*(i-1)+j; q{i1}=n2;
     n3=n3+mkron(q{:});
  end
  end

% initialize rising operators assuming level order: 1up,1down,2up,2down,...
%  * Sp for total spin SU(2)
%  * Qp for total charge SU(2)
%  * Cp for channel SU(3)
  clear Sp Qp Cp CP

  for i=1:NC
  for j=1:NC
   % acting on different spins (up and down!)
     i1=2*i-1; i2=2*j;
     q=repmat({e2},1,2*NC); q(min(i1,i2)+1:max(i1,i2)-1)={z2};
     q{i1}=f2';
     q{i2}=f2 ; Sp{i,j}=mkron(q{:});
     q{i2}=f2'; Qp{i,j}=mkron(q{:});

   % acting on the same spin (sum over spin later; for channel SU3)
     i1=2*i-1; i2=2*j-1;
     q=repmat({e2},1,2*NC); q(min(i1,i2)+1:max(i1,i2)-1)={z2};
     q{i1}=f2'; if i1~=i2
     q{i2}=f2 ; else q{i2}=f2'*f2; end
     Cp{i,j,1}=mkron(q{:});

     i1=2*i; i2=2*j;
     q=repmat({e2},1,2*NC); q(min(i1,i2)+1:max(i1,i2)-1)={z2};
     q{i1}=f2'; if i1~=i2
     q{i2}=f2 ; else q{i2}=f2'*f2; end
     
     Cp{i,j,2}=mkron(q{:});

     CP{i,j}=Cp{i,j,1} + Cp{i,j,2};
     if i==j
        if i==1 && j==1, d=size(Cp{1},1); end
      % make traceless
        CP{i,j}=CP{i,j}-diag(repmat(trace(CP{i,j})/d,1,d));
     end

  end
  end

% -------------------------------------------------------------------- %
% setup of symmetries
% -------------------------------------------------------------------- %
  clear IS q

  q.info='particle U(1)';
  q.type='A';
  q.Sp={};
  q.Sz={ sum(cat(3,CP{1:NC+1:end}),3) };
  q.qfac=[];
  q.jmap=[];
  IS=q;

  q.info='spin Sz U(1)';
  q.type='A';
  q.Sp={}; SP=sum(cat(3,Sp{1:NC+1:end}),3); % total Sp (sum over diagonal)
  q.Sz={ 0.5*comm(SP,SP') };
  q.qfac=2;
  IS(end+1)=q;

  q.info='channel SU(3)';
  q.type='SU3';
% NB! CP_13 = comm(CP_12,CP_23) already contained
  q.Sp={ CP{1,2}, CP{2,3} }; q.Sz={};
  for i=1:length(q.Sp), q.Sz{i}=0.5*comm(q.Sp{i},q.Sp{i}'); end
  q.Sz={ q.Sz{1}, q.Sz{1}+2*q.Sz{2} };
  q.qfac=2;               % factor(s) applied on z-labels
  q.jmap=[1 -0.5; 0 0.5]; % subsequent map for J-labels
  IS(end+1)=q;            % (multiplied from the right: J_out * qfac

% get proper symmetry eigenstates
  [U,Is]=getSymStates(IS);

  Is.qq=cat(1,Is.RR.J);
  QQ=uniquerows(chopd(Is.qq));
  Is.qz=mat2cell(Is.QZ,Is.dd,Is.d2);

  if size(QQ,1)~=size(Is.qq,1)
       wblog('WRN','some symmetry sectors not unique'); 
  else wblog(' * ','all symmetry sectors unique'); end

  fmt={}; inl(1); nI=numel(IS);

  for i=1:nI, n=numel(IS(i).Sz); l=5*n; m=length(IS(i).info);
     fprintf(1,sprintf(' %%%gs  ',l),IS(i).info);
     fmt{end+1}=repmat(' %4g',1,n);
     if l<m+3
     fmt{end}=[ fmt{end}, repmat(' ',1,m-l+3) ]; end
  end
  fmt=[ cat(2,fmt{:}) '\n' ]; inl(1);
  fprintf(1,fmt,QQ'); inl(1);

% transform ALL operators into new basis
  for i=1:numel(FF), FF{i}=U'*FF{i}*U; end
  for i=1:nI
     for j=1:numel(IS(i).Sz), IS(i).Sz{j}=U'*IS(i).Sz{j}*U; end
     for j=1:numel(IS(i).Sp), IS(i).Sp{j}=U'*IS(i).Sp{j}*U; end
  end

  SP=U'*SP*U;

% get F(2)-operator
  [F3u,Iu]=getSymmetryOps(FF{1,1},IS,FF(:,1));
  [F3d,Id]=getSymmetryOps(FF{1,2},IS,FF(:,2));

  F2=[QSpace(compactQS('A,A,SU3',Is.QZ,Is.QZ,Iu.QZ, cat(3,F3u{:})))
      QSpace(compactQS('A,A,SU3',Is.QZ,Is.QZ,Id.QZ, cat(3,F3d{:}))) ];

  N0=QSpace(contractQS(F2(1),[1 3],F2(1),[1 3])) + ...
            contractQS(F2(2),[1 3],F2(2),[1 3]);

  for i=1:2
     q=(1/3)*(QSpace(...
         contractQS(F2(i),[1 3],F2(i),[1 3])) + ...
         contractQS(F2(i),[2 3],F2(i),[2 3]));
     if ~isIdentityQS(q), error('Wb:ERR','invalid F2 operators!'); end
  end
  wblog('ok.','got valid F2 operator');

  E0=QSpace(getIdentityQS(F2));

  Z=N0+0*E0;
  for i=1:length(Z.data)
     z=round(diag(Z.data{i})); e=norm(diag(z)-Z.data{i});
     if e>1E-12, error('Wb:ERR','z-op not diagonal'); end
     Z.data{i}=diag((-1).^z);
  end

% spin operators (NB! spin behaves like J=1 spinor operator)
  sz=IS(2).Sz{1}; s2=0.5*(SP*SP' + SP'*SP)+sz*sz;
  S2=QSpace(compactQS('A,A,SU3',Is.QZ,Is.QZ,zeros(1,size(Is.QZ,2)), s2));
  Sz=QSpace(compactQS('A,A,SU3',Is.QZ,Is.QZ,zeros(1,size(Is.QZ,2)), sz));

  H0=-JH*S2 + B*Sz;

% operators for self-energy: FN = comm(H0,FC)
  FN=F2;
  for i=1:numel(FN)
      FN(i)=QSpace(contractQS(H0,2,FN(i),1))-contractQS(FN(i),2,H0,1,[1 3 2]);
      FN(i).info.otype=F2(i).info.otype;
  end
  FN=(1/JH)*FN; %  eliminate JH dependence (since H0 \ propto JH)

  A1=QSpace(getIdentityQS(F2,F2));   %   20 x 20 x    672 =>     64 x 64 x       4096
% A2=QSpace(getIdentityQS(A1,3,F2)); %  672 x 20 x  27456 =>   4096 x 64 x     262144
% A3=QSpace(getIdentityQS(A2,3,F2)); % 90?? x 10 x 3667?? => 262144 x 64 x 16.777.216

% generate auxiliary operator with empty state to generate A0
  Q=H0; Q.Q=repmat({zeros(1,size(A1.Q{1},2))},1,2);
  Q.data={1}; Q.info.cgs=repmat({1},1,3);
  A0=QSpace(getIdentityQS(Q,F2));

% ==================================================================== %
% NRG run
% ==================================================================== %

  if ~exist('Inrg','var'), cto lma ; param % kinit

     fout='./NRG/NRG';
     onrg={'fout',fout,'Nkeep',Nkeep}; inl(1);

     if exist('mat','var'), load2(mat); param
     else
      % having U(1)*U(1)*SU(3) => SU3
        mat=sprintf('%s_NRG_SU3_B%02.0f_D%d_L%02.0f_J%04.0f',...
        wbtstamp,-log10(B),Nkeep,Lambda*10,JH*1E4); param.mat=mat;

        if savepara, clear mat
        else
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           Inrg=NRGWilsonQS(H0,A0,Lambda,ff,F2,Z,onrg{:});
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           Imain.finished=datestr(now); clear NRGWilsonQS
           save([ mat '.mat'])

           EE=Inrg.EE; US=U; ES=E0; E0=Inrg.E0; clear U
           nrg_plot ; if isbatch, mfig([ mat '_nrg.pdf']), end
        end % of ~savepara
     end
  end

% op1=[];op2=F2; zflags=1;
% op1=[];op2=Sz; zflags=0;
% op1=[];op2=[F2;Sz]; zflags=[1 1 0];
% op1=[F2;FN;Sz]; op2=[F2;F2;Sz]; zflags=[1 1 1 1 0];
% op1=[F2(1);FN(1);Sz]; op2=[F2(1);F2(1);Sz]; zflags=[1 1 0];
  op1=[F2(1);FN(1);Sz]; op2=[F2(1);F2(1);Sz]; zflags=[1 1 0];

  o={'nostore','locRho'}; % ,'T',1E-4
  setopts(o,zflags,'T?');

  if savepara
     FC=F2; nostore=1; locRho=1;

     vars={sprintf('%s_%dops_setup.mat',param.mat(10:end),numel(op2)),...
       'param','H0','A0','Lambda','ff','fout','FC','FN','Z',...
       'op1','op2','Nkeep','nostore','locRho'};
     for v={'T','cflags','zflags'}
     if exist(v{1},'var'), vars{end+1}=v{1}; end, end

     f=[getenv('MEX') '/' vars{1}];
     fprintf(1,'\n   saving NRG environment to mat file\n');
     fprintf(1,'   mat: %s\n\n',repHome(f));
     save(f,vars{2:end});

     clear FC ; return
  end % savepara

  [om,a0,Idma]=fdmNRG_QS(fout,op1,op2, Inrg.Z,o{:});

  Imain.finished=datestr(now);
  save([ mat '.mat'])

  dma_plot ; if isbatch, mfig([ mat '_dma.pdf']), end

