
if 0
  cto lma
% load2 NRG_SU123/Wb110930_D10000_Etr70_l40_J0300 % th-ws-i7s06//LMA
  load2 NRG_SU123/Wb120117_D10000_Etr60_l40_j30  %th-ws-i710//LMA
else % TMPX

% clear ; savepara=1 ; rnrg_SUN
  setdef('savepara',0);

  Imain.cputime=cputime;

  inl(1);
  if ~exist('param','var')
     setdef('NC',3); Gamma=0.01;

   % Hund's coupling: 0.008 * 2.5/(3+0.5) * 4 = 0.0229
   % see rkondo.m for NC=3;
   % NB! factor 4 since Sop is defined without factor of 0.5 there!

   % JH=0.0023;
   % JH=0.0012; % too small (looses significance compared to Gamma!!)

   % JH=0.0300; % TK = ~3E-8 (Lambda=4 : converged at N>30; TK ~ k=25)
   % JH=0.0500; % TK = ~3E-12 (Lambda=4 : converged at N>45; TK ~ k=38)

 if 0
   % JH and Gamma *quite small* below => rescale units by factor 5
   % mixing small energies (JH,Gamma) with large band width
   % => truncation sets in too soon // Wb,Feb01,11
     Gamma=0.1; JH=0.22857142857143; % NB! S^2 op has another factor 1/4
 else
   % NB! rkondo => setupKondo_AM: 
   % H0 = J*SS with S in defined *without* factor 1/2)
     Gamma=0.01; JH=0.022857142857143;
   % =-4J = same as 'rkondo -k 3 4500' (abelian only)
 end

   % TST: Wb,Sep29,11
     Gamma=10; JH=30; emax=1E3; % LARGE_U
     ETRUNC=[repmat(max(10,10*JH),1,2),max(10,2*JH)];
     % ETRUNC(1): irrelevant
     % ETRUNC(2): skip high-lying states one iteration later ...
     % ETRUNC(3): clean cut within gap to higher-lying states

   % Nkeep=2048;
   % Lambda=3; L=40;
   % Nkeep=256; % TK factor of 10 too low (Nkeep too small)
   % Nkeep=260; % NB! 260 is min(Nkeep) to fully keep 2nd iteration
   % Nkeep=512;

   % Lambda=2.5; L=40;
   % NB! 260 is min(Nkeep) to fully keep 2nd iteration
   % Nkeep=260; % ac somewhat too high (ax itself fits pretty well, though!)
   % Nkeep=265; % ac slightly too high
   % Nkeep=266; % (same as 265)
   % Nkeep=267; % (same as 270)
   % Nkeep=268; % (same as 270)
   % Nkeep=270; % ac slightly too low
   % Nkeep=300; % ac somewhat too low

   % --------------------------------------------------------------------- %
   % rerun old rkondo NC=3 data % Wb,May29,08 => Wb,Sep27,11 (th-ws-i7s01)
   % TT=logspace(-8,-4,10); tfac=1;
   % Lambda=3; Etrunc=5; Nkeep=4096; L=40; seq=1;
   % --------------------------------------------------------------------- %
     TT=logspace(-8,-4,10); TT(1)=-1; tfac=1;

   % Nkeep=32; L=12; % for testing purposes
   % Nkeep=260; TT=logspace(-8,-4,128) % for testing purposes

   % ============================= SU123 ========================== runtime  IO-time %
   % Nkeep= 256; % NB! keeps up to  5070 states including z-labels    8 min
   % Nkeep= 512; % NB! keeps up to 10050 states including z-labels   22 min  >27 min!?
     % Nkeep=512 => 3.1G HD space @ 10Mb/sec => expecting 5 min!
     % (compare previous Abelian calculation: Nk=4096: runtime 20 hrs
     % eg. see $SCR/Data/goedel rkondo_070919_sge77429.log

   % Lambda=3; Etrunc=5; too small for given Lambda!
   % Lambda=4; Etrunc=5; ET1=-1; Nkeep=1E4; L=40; seq=1; % good
   % Lambda=4; Etrunc=6; Nkeep=4096; % ok, Nkeep=6000 is better

   % Lambda=4; Etrunc=6; ET1=-1; Nkeep=1E4; L=42; seq=1;
     % works nicely with LARGE_U [Nkept<1000!]
     % rho_disc ~ 1E-10 [ not quite converged yet ]
     % => need 1E-12 in agreement with Wb11_rho paper!
     Lambda=4; setdef('Etrunc',6); Nkeep=8000; L=42; seq=1;
     % works nicely with LARGE_U [Nkept<2100!]

   % Nkeep=6000; Lambda=4; L=40;   Etrunc=5.5;
   % Nkeep=6000; Lambda=3; L=50;   Etrunc=5.0;

   % Nkeep=2048; Lambda=4; L=40;   NKEEP=repmat(9086,1,4);
   % Nkeep=2048; Lambda=4; L=40;
   % Nkeep=1024; Lambda=4; L=40;
   % Nkeep= 512; Lambda=4; L=40;

   if isset('TSTflag')
     Nkeep= 260; Lambda=4; L=40; seq=1; clear Etrunc % nrgonly=1; % for testing purposes
   % Nkeep=  64; Lambda=4; L=20; seq=1; clear Etrunc % nrgonly=1; % for testing purposes
   end

   % Nkeep=2048; Lambda=3; L=50;
   % Nkeep=2560; Lambda=3; L=50;
   % Nkeep=3072; Lambda=3; L=50;
   % Nkeep=2560; Lambda=2.5; L=30;
   % Nkeep=3072; Lambda=2.5; L=50;

   % ============================= SU123 ========== for comparison: (SU2^4) %
   % Nkeep=1024; % NB! keeps up to 23500 states including z-labels
   % Nkeep=1500; % NB! keeps up to 31960 states including z-labels  (23500)
   % Nkeep=2048; % NB! keeps up to 53600 states including z-labels  (32800)
   % Nkeep=2600; % NB! keeps up to 63500 states including z-labels  (40400)
   % Nkeep=3500; % NB! keeps up to 63500 states including z-labels
   % ====================================================================== %

  end

  setupKondoAM % outsourced // Wb,Apr28,13

if isset('init_only'), inl(1); return, end
% if ~isbatch, wblog('==>','setup finished.'); inl(1); return, end

% ==================================================================== %
% NRG run
% ==================================================================== %

  if ~exist('Inrg','var'), cto lma ; param % kinit

     if exist('mat','var') && exist(mat,'file')
      % mat, exist(mat,'file'), which(mat)
        load2(mat); param
     else
      % setdef('fout','./NRG/NRG');

        tag=sprintf('D%04d',Nkeep); if exist('Etrunc','var') && Etrunc>0
        tag=[tag,sprintf('_Etr%2.0f',10*Etrunc)]; end
        tag=[tag,sprintf('_l%02.0f_j%02.0f',10*Lambda,10*JH/Gamma)];

        mat='NRG_SU123';
        fout=sprintf('%s/%s/NRG',mat,tag);
        mat=sprintf('%s/%s_%s',mat,wbtstamp,tag);

        [i,istr]=system(['mkdir.pl ' fout]); if i, error('Wb:ERR','???'); end

      % onrg={'fout',fout,'Nkeep',Nkeep};
        onrg=setopts('-',fout,Nkeep,'NKEEP?','ETRUNC?','Etrunc?','ET1?');
        inl(1); dispopts(onrg);

        if NC==1, setopts(onrg,'zflag',2); end

      % having U(1)*SU(2)*SU(3) => SU123
        if savepara, clear mat
        else
         % ============================================= %
           Inrg=NRGWilsonQS(H0,A0,Lambda,ff,FF,Z,onrg{:});
           Imain.finished=datestr(now); clear NRGWilsonQS
           save([ mat '.mat'])
         % ============================================= %

           EE=Inrg.EE; E0=Inrg.E0;
           nrg_plot ; if isbatch, mfig([ mat '_nrg.pdf']), end

           II=getrhoNRG_resall(fout,'-v');
           if isbatch, mfig([mat '_rhores.pdf']); end
           save([ mat '.mat'],'II','-append')

        end % of ~savepara
     end
  end

  if exist('nrgonly','var') && nrgonly, return; end

end % of TMPX

% banner('no nostore!'); odma={}; op2=FF; zflags=1;
  odma={'nostore','locRho'}; setopts(odma,'T?','emax?');

  if savepara
     FC=FF; nostore=1; locRho=1;

     vars={sprintf('%s_setup.mat',param.mat),...
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

  if ~exist('tfac','var'), tfac=0.5; end

% used by dma_plot
  yli=[-0.1 1.2];

if exist('TT','var'), P=PSet('T',TT); else P=PSet('T',-1); end
for ip=1:P.n, [p,pstr,tstr]=P(ip); structexp(p); 
  banner('%s\n%s',pstr,tstr);

% T<0 takes default in fdmNRG_QS()!
  if T>0, eps=tfac*T; else eps=[]; end
  clear om a0 Idma

  if exist('nchi','var') && ~isempty(nchi) && ip>nchi
     op1_=op1; op1=op1(1:end-1);
     op2_=op2; op2=op2(1:end-1);
  end

  if seq
     if isempty(op1), op1=op2; end
     for io=1:numel(op2)
        [om{io},a0{io},Idma(io)]=fdmNRG_QS(fout,op1(io),op2(io), ...
        Inrg.ops.Z,odma{:},'zflags',zflags(io),'T',T);
     end
     if norm(diff(cat(2,om{:}),[],2))>1E-12
        error('Wb:ERR','omega mismatch'); end
     om=om{1}; a0=cat(2,a0{:});
  else
     setopts(odma,zflags);
     [om,a0,Idma]=fdmNRG_QS(fout,op1,op2,Inrg.ops.Z, odma{:});
  end

  if T<0, TT(ip)=Idma.T; end

  TK = 1./(4*sum(Idma(end).reA0(end-1:end)));
  wblog('TST','TK=%.3g',TK);

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

% rerun getrhoNRG_res to blend in data from dma_plot
  II=getrhoNRG_resall(fout,'-v');
  if isbatch, mfig([mat '_rhores.pdf'],'-f'); end
  save([ mat '.mat']);

