
% test routine for NRGWilsonQS.cc using symmetries
% Wb,Sep08,17

  if ~exist('U','var') || ~exist('param','var') || isempty(param)
     U=0.12; epsd=-U/2; Gamma=0.01; B=0; Lambda=2;
     N=55; Nkeep=512;
  end

  global param
  param=add2struct('-',U,epsd,Gamma,B,N,Lambda);

  [FF,Z,SS,IS]=getLocalSpace('FermionS','Acharge,Aspin','-v');

  N0=FF;
     for i=1:numel(FF)
     N0(i)=contractQS(FF(i),'13*',FF(i),'13'); end
  N0=sum(N0);

  H0 = epsd*N0 + B*SS(3) + (U/2)*N0*(N0-IS.E) + 0*IS.E;
  A0 = QSpace(permuteQS(getIdentityQS(getvac(IS.E),IS.E),[1 3 2]));

  ff=getNRGcoupling(Gamma, Lambda, N);
  ff=ff.* (Lambda.^((1:(N-1))/2));

  HK=H0; AK=A0; FKK=QSpace(1,2); rL=sqrt(Lambda);
  fprintf(1,'\n'); clear E0; nn=nan(N,2); EX=nan(N,2); At=QSpace;

  setdef('sflag',1);
  if sflag
     fout='NRG/NRGX'; cto lma
     wblog('I/O','using %s_##.mat',fout); 
  end

  structdisp(param)

  for k=1:N
     if k<N || ~sflag
           o={'Nkeep',Nkeep};
     else, o={'Nkeep',0}; end

     [ee,I]=eigQS(HK,o{:}); ee=ee(:,1);

     if ~isempty(I.EK.data)
        q=getDimQS(I.EK); nn(k,1)=q(2);
        EX(k,1)=max([I.EK.data{:}]);
     end
     if ~isempty(I.ED.data)
        q=getDimQS(I.ED); nn(k,2)=q(2);
        EX(k,2)=min([I.ED.data{:}]);
     end

     if k==1, E0=0; else E0(k)=min(ee); end;
     ee=ee-E0(k);

     if sflag
     At=contractQS(AK,2,I.AD,1,[1 3 2]); end
     Ak=contractQS(AK,2,I.AK,1,[1 3 2]);

     Ek=QSpace(I.EK); if ~isempty(Ek.data), Ek=Ek-E0(k); end
     Et=QSpace(I.ED); if ~isempty(Et.data), Et=Et-E0(k); end
     if sflag
        q=struct('AK',Ak,'AD',At,'HK',Ek,'HD',Et,'E0',E0(k));
        save(sprintf('%s_%02g.mat',fout,k-1),'-struct','q');
     end

     m=min(Nkeep,length(ee));
     EE(1:m,k)=ee(1:m);

     if k==N, break; end

     for s=1:2
        FKK(s)=contractQS(Ak,'13*',contractQS(Ak,3,FF(s),2),'13');
     end

     AK=permute(QSpace(getIdentityQS(I.AK,2,IS.E,2)),[1 3 2]);

     HK=diag(Ek*rL);
     HK=QSpace(contractQS(AK,'13*',contractQS(HK,2,AK,1),'13'));

     for s=1:2
        Q=QSpace(contractQS(AK,'13*',contractQS(FKK(s),'13*', ... % = FKK'
          contractQS(AK,3,ff(k)*Z*FF(s),2),'14'),'13'));
        HK=HK+Q+Q';
     end

     fprintf(1,' %4d/%d (%g) EK=%.3g ... \r',k,N,...
       max(getDimQS(HK)),EX(k,1));

     HK=HK+1E-99*QSpace(getIdentityQS(AK,2));
  end

  fprintf(1,'\n\n');

  Inrg.E0=E0;
  Inrg.EScale=Lambda.^(-(0:length(E0)-1)/2);
  Inrg.phE0=sum(E0.*Inrg.EScale);
  Inrg.NK=nn;
  Inrg.N=N;

  EE(find(EE==0))=nan; EE(1,find(isnan(EE(1,:))))=0;

  if ~exist('plotflag','var') || plotflag
     param.D=Nkeep;
     param.L=N;
     nrg_plot; set(ah(1:2),'YLim',[0 6]);
     drawnow
  end

if ~sflag, return, end

  Inrg.Lambda=Lambda;
  f=sprintf('%s_info.mat',fout);
  save(f,'-struct','Inrg');

return

  op1=[]; op2=FF(1); zflags=1;
  nostore=1; locRho=1;

  odma=setopts('-','T?','zflags?','cflags?','rhoNorm?',...
   'nlog?','emin?','emax?','-calcOps?','-nostore?','-locRho?');

  [om,a0,Idma] = fdmNRG_QS(fout,op1,op2,Z, odma{:});
  dma_plot

