
% setenv CG_VERBOSE 6

if 0
  cto lma
% load2 NRG_Sp23/Wb120119_D8000_Etr70_l40_j30.mat  % th-ws-i710//LMA
% load2 NRG_Sp23/Wb120123_D4096_Etr70_l40_j30.mat  % th-ws-i710//LMA
% load2 NRG_Sp23/Wb120124_D2048_Etr50_l40_j30.mat  % th-ws-i710//LMA
% load2 NRG_Sp23/Wb120212_D8000_Etr70_l20_j30.mat  % th-cl-i7tb01//LMA
  load2 NRG_Sp23/Wb120523_D4000_Etr70_l40_j30.mat  % th-ws-i7s08//LMA
  mat=[mat '-2']; seq=0; clear nostore
else % RUN_NRG

% clear ; savepara=1 ; rnrg_SpN
  setdef('savepara',0);

  Imain.cputime=cputime;

  inl(1);
  if ~exist('param','var')
     NC=3; Gamma=0.01;

   % Gamma=0.01; JH=0.022857142857143;

   % TST: Wb,Sep29,11
     Gamma=10; JH=30; emax=1E3; % LARGE_U
     ETRUNC=[repmat(max(10,10*JH),1,2),max(10,2*JH)];
     % ETRUNC(1): irrelevant
     % ETRUNC(2): skip high-lying states one iteration later ...
     % ETRUNC(3): clean cut within gap to higher-lying states

   % expected TKondo (see tex notes)
     global g_TKondo
     g_TKondo=sqrt(JH*Gamma/2)*exp(-(NC+1)*JH/Gamma); % 4=2S+1 = 2(S+1/2)
     
   % Lambda=4; Etrunc=5; Nkeep=2048; L=42; seq=0;
     Lambda=4; setdef('Etrunc',6); Nkeep=4096; L=42; seq=0;
   % Lambda=2; Etrunc=7; Nkeep=8000; L=84; seq=1;

     TT=logspace(-8,-4,10); TT(1)=-1; tfac=1;

     global param
     param=add2struct('-',Gamma,JH,Lambda,L,NC,Nkeep,...
       'NKEEP?','Etrunc?','ETRUNC?');

     ff=getNRGcoupling(Gamma,Lambda,L);
     param.ff=ff;
  end

% wbc=wbtoc;
% start(wbc,'getLocalSpace');

% setting up SU(3) symmetric Anderson/Kondo hybrid
% Wb,Sep26,11
  [FF,Z,S3,IS]=getLocalSpace('FermionS','SU2spin,SpNchannel','NC',NC,'-v');
% stop(wbc,'getLocalSpace');

% start(wbc,'contractQS');
  S2=QSpace(contractQS(S3,[1 3],S3,[1 3],'conjA'));
% stop(wbc,'contractQS');
  H0=-JH*S2;

% operators for self-energy: FN = comm(H0,FC)
  FN=FF;
  
% resume(wbc,'contractQS');
  FN=QSpace(contractQS(H0,2,FN,1))-contractQS(FN,2,H0,1,[1 3 2]);
% stop(wbc,'contractQS');
  FN.info.otype=FF.info.otype;
  FN=(1/JH)*FN; %  eliminate JH dependence (since H0 \ propto JH)

% start(wbc,'getIdentityQS');
  A1=QSpace(getIdentityQS(FF,FF)); if exist('init_only','var') && init_only
  A2=QSpace(getIdentityQS(A1,3,FF)); %   61 x 4 x   1232 =>   4096 x 64 x    262.144
% A3=QSpace(getIdentityQS(A2,3,FF)); % 1232 x 4 x  31640 => 262144 x 64 x 16.777.216
% A2 takes about 28 seconds; generates up to d=614 dimensional multiplets (7MB)
% A3 takes about 1/2 hr; generates up to d=2240 dimensional multiplets (340MB)
% stop(wbc,'getIdentityQS');
  end

% resume(wbc,'getIdentityQS');
  A0=permute(QSpace(getIdentityQS(getvac(Z),FF)),[1 3 2]);
% stop(wbc,'getIdentityQS');

% wbc

if isset('init_only'), inl(1); return, end
if ~isbatch, wblog('==>','setup finished.'); inl(1); return, end

% ==================================================================== %
% NRG run
% ==================================================================== %

  if ~exist('Inrg','var'), cto lma ; inl(1); disp(param) % kinit

     if exist('mat','var') && exist(mat,'file')
      % mat, exist(mat,'file'), which(mat)
        load2(mat); param
     else
      % setdef('fout','./NRG/NRG');

        tag=sprintf('D%04d',Nkeep); if exist('Etrunc','var') && Etrunc>0
        tag=[tag,sprintf('_Etr%2.0f',10*Etrunc)]; end
        tag=[tag,sprintf('_l%02.0f_j%02.0f',10*Lambda,10*JH/Gamma)];

        mat='NRG_Sp23';
        fout=sprintf('%s/%s/NRG',mat,tag);
        mat=sprintf('%s/%s_%s',mat,wbtstamp,tag);

        [i,istr]=system(['mkdir.pl ' fout]); if i, error('Wb:ERR','???'); end

      % onrg={'fout',fout,'Nkeep',Nkeep};
        onrg=setopts('-',fout,Nkeep,'NKEEP?','ETRUNC?','Etrunc?','ET1?');
        dispopts(onrg);

      % having U(1)*SU(2)*SU(3) => SU123
        if savepara, clear mat
        else
         % ======================================================== %
           Inrg=NRGWilsonQS(H0,A0,Lambda,ff,FF,Z,onrg{:},'zflag',2);
           Imain.finished=datestr(now); clear NRGWilsonQS
           save([ mat '.mat'])
         % ======================================================== %

           EE=Inrg.EE; E0=Inrg.E0;
           nrg_plot ; if isbatch, mfig([ mat '_nrg.pdf']), end

           II=getrhoNRG_resall(fout,'-v');
           if isbatch, mfig([mat '_rhores.pdf']); end
           save([ mat '.mat'],'II','-append')

        end % of ~savepara
     end
  end

  if exist('nrgonly','var') && nrgonly, return; end

end % of RUN_NRG

  if exist('NKEEP','var'), n=max([NKEEP Nkeep]); else n=Nkeep; end
  if size(Inrg.NK,2)>=4, n=max(Inrg.NK(:,1)); end

  setdef('seq',0);
  if n<2048 || seq
     wblog(' * ','using 3 operators (n=%g)',n); 
   % op1=[]; op2=[FF,S3]; zflags=[1 0]; % kinit
     op1=[FF;FN;S3];
     op2=[FF;FF;S3]; zflags=[1 1 0];
  else
     op1=[]; op2=FF; zflags=1;
  end

  nostore=0;

% banner('no nostore!'); odma={}; op2=FF; zflags=1;
  odma={'locRho'}; setopts(odma,'T?','emax?','-nostore?');
  setdef('tfac',0.5);

% used by dma_plot
  yli=[-0.1 1.2];

if exist('TT','var'), P=PSet('T',TT); else P=PSet('T',-1); end

try
for ip=1:P.n, [p,pstr,tstr]=P(ip); structexp(p); 
  banner('%s\n%s',pstr,tstr);

% T<0 takes default in fdmNRG_QS()!
  if T>0, eps=tfac*T; else eps=[]; end
  clear om a0 Idma

  if exist('nchi','var') && ~isempty(nchi) && ip>nchi
     op1_=op1; op1=op1(1:end-1);
     op2_=op2; op2=op2(1:end-1);
  end

  if isset('seq')
     if isempty(op1), op1=op2; end
     for io=1:numel(op2)
        [om{io},a0{io},Idma(io)]=fdmNRG_QS(fout,op1(io),op2(io), ...
        Inrg.ops.Z,odma{:},'zflags',zflags(io),'T',T);

        save([ mat '.mat']); % clear function fdmNRG_QS
     end

     if norm(diff(cat(2,om{:}),[],2))>1E-12
        error('Wb:ERR','omega mismatch'); end
     om=om{1}; a0=cat(2,a0{:});
  else
     setopts(odma,T,zflags);
     [om,a0,Idma]=fdmNRG_QS(fout,op1,op2,Inrg.ops.Z, odma{:});

     save([ mat '.mat']); % clear function fdmNRG_QS
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

catch
  wblog('WRN',''); whos; l=lasterror, disp(l.message), dispstack(l.stack);
  save([ mat '-v7.3.mat'],'-v7.3'); wblog('TST',''); clear fdmNRG_QS
  save([ mat '.mat']); wblog('WRN','');
end

% rerun getrhoNRG_res to blend in data from dma_plot
  II=getrhoNRG_resall(fout,'-v');
  if isbatch, mfig([mat '_rhores.pdf'],'-f'); end
  save([ mat '.mat']);

