% SIAM impurity parameters (all energies in units of the half bandwidth D:=1)
%
%   U        onsite interaction
%   epsd     impurity level position
%   Gamma    hybridization strength of impurity
%   B        magnetic field acting on impurity ($-B S_z$; Zeeman splitting)

  setdef('U',0.12,'Gamma',0.01,'B',0);
  setdef('epsd',-U/2);

% NRG discretization parameters
%
%   Lambda   coarse graining strength
%   N        length of Wilson chain
%   Nkeep    number of multiplets (or states if Abelian) to keep
%   Etrunc   energy truncation threshold (complementary to Nkeep)
%   z        shift of logarithmic discretization, as in $\Lambda^{(z-n)/2}$ with $z\in[0,1[$
%
% By setting default Nkeep large, preference is on truncation by energy;
% using Nkeep nevertheless, thus serves as safeguard

  setdef('Lambda',2,'N',55,'Nkeep',1024,'Etrunc',7,'z',0);

  global param
  param=add2struct('-',U,epsd,Gamma,B,N,Lambda,z);

  ff=getNRGcoupling(Gamma,Lambda,N,'z',z);

  if ~isset('B'), SYM='Acharge,SU2spin'; else SYM='Acharge,Aspin'; end
  [FF,Z,SS,IS]=getLocalSpace('FermionS',SYM,'-v');

  N0=FF; nF=numel(FF);
     for i=1:nF
     N0(i)=contract(FF(i),'13*',FF(i),'13'); end
  N0=sum(N0);

  A0 = getIdentity(getvac(IS.E),IS.E,[1 3 2]);

  H0 = epsd*N0 + (U/2)*N0*(N0-IS.E) + 0*IS.E;
  if isset('B'), H0 = H0 - B*SS(end); end

  rL=sqrt(Lambda);  % rL = `root' Lambda
  ff=ff.* rL.^(1:N-1);

  setdef('sflag',1);
  if sflag
     fout={'NRG' '/myNRG'}; cto lma; 
     if ~isdir(['./' fout{1}]), system(['mkdir ' fout{1}]); end
     fout=[fout{:}];
     wblog('I/O','using %s_##.mat',fout); 
  end

% initialize iterative diagonalization
%    K  kept space      (AK $\equiv\ A_K$)
%    D  discarded space (AD $\equiv\ A_D$)
% scalar operators only have KK $\to$ K, DD $\to$ D (like HK $\equiv\ H_{K}$, HD $\equiv\ H_{D}$)
% general operators also have off-diagonal blocks KD, DK (cf. FKK $\equiv\ F_{KK}$)
  AD=QSpace; k=0;
  AK=setitags(setitags(A0,'-A',k),1,'Lvac');    % sets itags {'Lvac','K00','s00'}
  HK=setitags(H0,'-op:K',k);                    % sets itags {'K00','K00*'}

  FKK=QSpace(1,nF);

  EX=nan(N,2);
  nn=nan(N,2);
  E0=zeros(1,N);

  Inrg=struct('N',N,'Lambda',Lambda,'EK',QSpace(1,N));

  fprintf(1,'\n'); disp(param);
  fprintf(1,'\n');

  for k=1:N
     if k<N || ~sflag
          o={'Nkeep',Nkeep,'Etrunc',Etrunc};
     else o={'Nkeep',0};
     end

     [ee,I]=eigQS(HK,o{:});
     ee=ee(:,1);
     EK=QSpace(I.EK);
     ED=QSpace(I.ED);

     if k>1 E0(k)=min(ee); end
     ee=ee-E0(k);

     if EK
        Dk=getDimQS(EK); nn(k,1)=Dk(1,2);
        EK=EK-E0(k); q=EK.data;
        EX(k,1)=max([q{:}]);
     end
     if ED
        q=getDimQS(ED); nn(k,2)=q(1,2);
        ED=ED-E0(k); q=ED.data;
        EX(k,2)=min([q{:}]);
     end

     if sflag
     AD=contract(AK,I.AD,[1 3 2]); end
     AK=contract(AK,I.AK,[1 3 2]);
     if sflag
        q=struct('AK',AK,'AD',AD,'HK',EK,'HD',ED,'E0',E0(k));
        save(sprintf('%s_%02g.mat',fout,k-1),'-struct','q');
     end

     m=min(Nkeep,length(ee));
     EE(1:m,k)=ee(1:m);
     Inrg.HK(k)=EK;

     if k==N, break; end

     for i=1:nF
        FKK(i)=contract(AK,'!2*',{AK,FF(i),'-op:^s'});
     end

     AK=getIdentity(I.AK,2,IS.E,[1 3 2]);
     AK=setitags(AK,'-A',k);

     HK=diag(EK*rL);
     HK=contract(AK,'!2*',{HK,AK});

     for i=1:nF
        q=ff(k)*Z*FF(i);
        Q=contract(AK,'!2*',{FKK(i),'*',{AK,q,'-op:^s'}});
        HK=HK+Q+Q'; % +Q' adds hermitian conjugate (" + H.c.")
     end

     fprintf(1,' %4d/%d (%g) EK=%.3g ... \r',k,N,max(Dk(1,:)),EX(k,1));

     HK=HK+0*getIdentity(AK,2);

  end
  fprintf(1,'\n\n');

  Inrg.E0=E0;
  Inrg.EScale=Lambda.^(-(0:length(E0)-1)/2);
  Inrg.phE0=sum(E0.*Inrg.EScale);
  Inrg.NK=nn;

  if sflag
     f=sprintf('%s_info.mat',fout);
     save(f,'-struct','Inrg');
  end

  EE(find(EE==0))=nan;
  EE(1,find(isnan(EE(1,:))))=0;

  if ~exist('plotflag','var') || plotflag
     param.D=Nkeep;
     param.L=N;

     nrg_plot

     if isset('Etrunc') && Etrunc>4, q=Etrunc-1; else q=5; end
     set(ah(1:2),'YLim',[0 q]);

     set(ah(4),'YLim',[0 q+1]);
     set(ah(5),'YLim',[0 max(Inrg.NK(:))]);

     h=header('NW'); header('NE',''); header('SE','')
     set(h,'String',regexprep(get(h,'String'),'NRGWilsonQS','xnrgQS'));
     mfig
  end

  return

  if sflag
     op1=[];
     op2=FF(1);
     zflags=1;

     nostore=1;
     locRho=1;

     ofdm=setopts('-','T?','zflags?',...
      'nlog?','emin?','emax?','-nostore?','-locRho?');

     [om,a0,Ifdm] = fdmNRG_QS(fout,op1,op2,Z, ofdm{:});
     fdm_plot
  end

