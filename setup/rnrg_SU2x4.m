
% cto lma % Wb,Oct01,11
% load NRG_SU2x4/Wb110930_D9000_Etr60_l40_J300000 % middlehead
% TT=logspace(-8,-4,10); TT(1)=-1; tfac=1; seq=1;
% if ~exist('S3','var') || isempty(S3), S3=SS; clear ss; end

  Imain.cputime=cputime;

% -------------------------------------------------------------------- %
  if ~exist('param','var') || ~exist('Inrg','var'), inl(1);
% tag: RNRG
% -------------------------------------------------------------------- %

  setdef('NC',3); Gamma=0.01; 

% JH=0.023; % Hund's coupling
% JH=0.03;

  if 0
    % JH and Gamma *quite small* below => rescale units by factor 5
    % mixing small energies (JH,Gamma) with large band width
    % => truncation sets in too soon // Wb,Feb01,11
      Gamma=0.1; JH=0.22857142857143; % NB! S^2 op has another factor 1/4
  else
    % NB! rkondo => setupKondo_AM: 
    % H0 = J*S2 with S in defined *without* factor 1/2)
      Gamma=0.01; JH=0.022857142857143;
    % =-4J = same as 'rkondo -k 3 4500' (abelian only)
  end

    % TST: Wb,Sep29,11
      Gamma=10; JH=30; emax=1E3; % LARGE_U
      ETRUNC=[repmat(max(10,10*JH),1,2),max(10,2*JH)];
      % ETRUNC(1): irrelevant
      % ETRUNC(2): skip high-lying states one iteration later ...
      % ETRUNC(3): clean cut within gap to higher-lying states

  Lambda=4; L=40; seq=1;

% Nkeep=128; % for testing purposes
% Nkeep=388; % includes 3rd site exactly

% NB! equiv to keeping 23500 states including z-labels
% Nkeep=1500; % 23500
% Nkeep=2048; % 32800
% Nkeep=2600; % 43220
% Nkeep=3072; %      

% Nkeep=9000; Etrunc=6; % using LARGE_U
  Nkeep=9000; Etrunc=7; % Etrunc=6 leads to converged results for SU123

% Nkeep=2600; Etrunc=4;
% Nkeep=9000; Etrunc=5; ET1=-1;
% Nkeep=5000; Etrunc=6; # running large (unfinished FDM)
% Nkeep=5000; Etrunc=8; # running large (stopped already during NRG)

  TT=logspace(-8,-4,10); TT(1)=-1; tfac=1;

  global param
  param=add2struct('-',Gamma,JH,Lambda,L,NC,Nkeep,'NKEEP?','Etrunc?');

  ff=getNRGcoupling(Gamma,Lambda,L);
  param.ff=ff;

% setting up SU(2)^4 for symmetric Anderson/Kondo hybrid
% Wb,Jun10,10

% for last local setup, see Archive/rnrg_SU2x4_110831.m
  [FF,Z,S3,IS]=getLocalSpace('FermionS','SU2spin,SU2charge(:)','NC',NC,'-v');

  S2=QSpace(contractQS(S3,[1 3],S3,[1 3]));

% auxiliary identity operator for EMPTY state to generate A0
  Q=Z; Q.Q=repmat({zeros(1,size(Z.Q{1},2))},1,2); Q.data={1};
  Q.info.cgs=repmat({1},1,size(Q.info.cgs,2));

  A0=QSpace(getIdentityQS(Q, 1,Z,1));
  A1=QSpace(getIdentityQS(A0,3,Z,1));
  A2=QSpace(getIdentityQS(A1,3,Z,1)); % still fast

% NB! A3 requires 24G (!!) leads to swapping!
% * A3 on its own is still doable using save2Mx()
%   instead of toMx() in getIdentityQS to avoid additional copy
% * 8214 blocks with total final dimensions
%   14229 x 13 x 590856 => 262144 x 64 x 16,777,216 = 1.7E+07
%   => average gain of factor of 28 in final matrix dimension
% A3=QSpace(getIdentityQS(A2,3,Z,1)); % takes about 40 secs

% build H0
  H0=-JH*S2;

% test first Wilson iteration
% NB! H still shows plenty of degeneracy (even so for H0!)
  H=QSpace(contractQS(A1,[1 2],contractQS(H0,2,A1,1),[1 2]));
  for i=1:numel(FF)
     Q=contractQS(FF(i),1,A1,1,'-L'); % NB! LsR symmetry of A1
     Q=contractQS(ff(1)*FF(i),[3 2],Q,[2 3]);
     H=H+contractQS(A1,[1 2],Q,[2 1]);
  end

% operators for self-energy: FN = comm(H0,FC)
  FN=FF;
  for i=1:numel(FF)
     FN(i)=QSpace(contractQS(H0,2,FF(i),1))-contractQS(FF(i),2,H0,1,[1 3 2]);
     FN(i).info.otype=FF(i).info.otype;
  end
  FN=(1/JH)*FN; %  eliminate JH dependence (since H0 \ propto JH)

if isset('init_only'), inl(1); return, end
if ~isbatch, wblog('==>','setup finished.'); inl(1); return, end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %
% NRG run

  Imain.started=datestr(now);
  cto lma ; param %  kinit

  if exist('mat','var') && exist(mat,'file')
   % mat, exist(mat,'file'), which(mat)
     load2(mat); param
  else
     tag=sprintf('D%04d',Nkeep); if exist('Etrunc','var') && Etrunc>0
     tag=[tag,sprintf('_Etr%2.0f',10*Etrunc)]; end
     tag=[tag, sprintf('_l%02.0f_j%02.0f',10*Lambda,10*JH/Gamma)];

     mat='NRG_SU2x4';
     fout=sprintf('%s/%s/NRG',mat,tag);
     mat=sprintf('%s/%s_%s',mat,wbtstamp,tag);

     [i,istr]=system(['mkdir.pl ' fout]); if i, error('Wb:ERR','???'); end
  end

  onrg=setopts('-',fout,Nkeep,'NKEEP?','Etrunc?','ETRUNC?','ET1?');
  inl(1);

% ********************************************************* %
  Inrg=NRGWilsonQS(H0,A0,Lambda,ff,FF,Z,onrg{:},'zflag',2);
% ********************************************************* %
  Imain.finished=datestr(now); clear NRGWilsonQS
  save([ mat '.mat'])

  EE=Inrg.EE; E0=Inrg.E0;
  nrg_plot ; if isbatch, mfig([ mat '_nrg.pdf']); end

  II=getrhoNRG_resall(fout,'-v');
  if isbatch, mfig([mat '_rhores.pdf']); end
  save([ mat '.mat'],'II','-append')

% -------------------------------------------------------------------- %
  end % of RNRG
% -------------------------------------------------------------------- %

  if exist('nrgonly','var') && nrgonly, return; end
% return

  if exist('NKEEP','var'), n=max([NKEEP Nkeep]); else n=Nkeep; end

  setdef('seq',0);
  if n<2100 || seq
   % op1=[]; op2=[FF,S3]; zflags=[1 0]; % kinit
     op1=[FF(1);FN(1);S3];
     op2=[FF(1);FF(1);S3]; zflags=[1 1 0];
  else
     op1=[]; op2=FF(1); zflags=1;
  end

  o={'nostore','locRho'}; setopts(o,'emax?');

if exist('TT','var'), P=PSet('T',TT); else P=PSet('T',-1); end
for ip=1:P.n, [p,pstr,tstr]=P(ip); structexp(p); 
  banner('%s\n%s',pstr,tstr);

% T<0 takes default in fdmNRG_QS()!
  if T>0, eps=tfac*T; else eps=[]; end
  clear om a0 Idma

  setopts(o,'T?','-q');

  if seq
     if isempty(op1), op1=op2; end
     for io=1:numel(op2)
        [om{io},a0{io},Idma(io)]=fdmNRG_QS(fout,op1(io),op2(io), ...
        Inrg.Z, o{:},'zflags',zflags(io));
     end
     if norm(diff(cat(2,om{:}),[],2))>1E-12
        error('Wb:ERR','omega mismatch'); end
     om=om{1}; a0=cat(2,a0{:});
  else
     setopts(o,zflags);
     [om,a0,Idma]=fdmNRG_QS(fout,op1,op2,Inrg.Z, o{:});
  end

  if T<0, TT(ip)=Idma.T; end

  TK = 1./(4*sum(Idma(end).reA0(end-1:end)));
  wblog('TST','TK=%.3g',TK);

  if size(a0,2)>=2 && abs(sum(a0(:,2)))<1E-6
   % also show improved spectral function (assuming that a0(:,2) is self-energy)
     ac_cmd='[ac,Ic]=getGC_JJH(ox,ax(:,1),-JH*ax(:,2),''-q'');';
  else clear ac_cmd; end

  IDMA(ip,1:numel(Idma))=Idma; Idma=Idma(1);
  dma_plot

  IALL(ip)=add2struct('-',om,a0,ox,ax,TK);

  Imain.finished=datestr(now);
  Imain.cputime=cputime - Imain.cputime;
  save([ mat '.mat'])

  if size(a0,2)==3 && abs(sum(a0(:,3)))<1E-6
     a0_=a0; a0=a0(:,1:2); % dma_plot without chi
     dma_plot
  end

  if isbatch || P.n>1
    if P.n>1, m=sprintf('%s_%02g_dma.pdf',mat,ip);
    else m=sprintf('%s_dma.pdf',mat); end
    mfig(m,'-f','-q');
  end

  sysinfo

end % of ip

% rerun gethroNRG_res to blend in data from dma_plot
  II=getrhoNRG_resall(fout,'-v');
  if isbatch, mfig([mat '_rhores.pdf'],'-f'); end
  save([ mat '.mat']);

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

