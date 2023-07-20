function [r2,E0,Il,HAM,RR,EE]=sweepPsi(HAM,varargin)
% function [r2,E0,Il,HAM,RR,EE]=sweepPsi(HAM [,opts])
%
%    Perform single DMRG sweep, starting from kc.
%
%    All provided options (except for the options -[vtT])
%    are directly handed over to update_psi_2site().
%    This includes the flag '-b' for a plain bond update,
%    instead of the default enlarged bond update (which
%    corresponds to a full 2-site update!).
%
% Options
%
%    -v    more verbose output
%
%    -T    compute ground state energy in single sweep
%          through the system *without* optimization
%          i.e. starting from kc, and returning to kc
%          e.g. to generate e0 data, etc.
%    -t    same as -T, but continue with actual DMRG sweep.
%
% Example usage
%
%    [r2,E0]=sweepPsi(HAM,'-T');
%
% AW (C) 2015

% added possibility to keep NPsi states within the same symmetry sector
% for former version, see: Archive/initNRG_160728.m
% Wb,Jul06,16 // Kyoto

% TST clear; wsys='SU2_w2'; NPsi=1; dbwrn; tst_Hamilton1D
% TST clear; wsys='SU2_w2'; NPsi=4; Qtot=[]; dbwrn; tstflag=0; L=10; tst_Hamilton1D
% TST clear; wsys='SU2_w2'; NPsi=8; Qtot=[]; dbwrn; rtol=1E-20; L=10; tst_Hamilton1D

  if nargin<1 || nargout>6
    helpthis, if nargin || nargout, wbdie('invalid usage'), end
    return
  end

  getopt('INIT',varargin);
     vflag=getopt('-v');

     if getopt('-t'), tflag=1;
     elseif getopt('-T'), tflag=2;
     else tflag=0; end

  o2site=getopt('get_remaining');

  if vflag, o2site{end+1}='-v'; end

  X1=load_dmrg_data(HAM,1);
  if isempty(X1.AK)

     if 1
        kc=initHK(HAM);
     else
        kc=initPsi(HAM);
        k2=initHK(HAM,kc);
        if k2~=kc, wbdie('inconsistent kc setting (%g/%g)',kc,k2);
        end
     end

  else kc=getCurrentSite(HAM); end

  L=numel(HAM.mpo);

  if tflag
     [E,e0,r2]=getEnergies(HAM,kc); E=E/L; ee=E(find(~isnan(E)));
     e=max(abs(diff(ee)));
     if e>1E-12, wbdie('got severe energy inconsistency (@ %.3g)',e); 
     else
        e=[mean(ee), std(ee)];
        wblog('ok.',['got consistent energy over MPS ' ...
       '(%.6g @ %.3g)'],e(1),e(2)/abs(e(1)));
     end

     if tflag>1, Il=e0; return, end
  end

  k1=1; k2=L; skipped=0;
  if kc~=1
     wblog('WRN','continuing at current site kc=%g',kc);
     k1=max(1,kc); skipped=k1-1; printf('\n')
  end

  s={'DMRG sweep'}; 

  Nkeep = getopt(HAM.info,'sweep','Nkeep',[]);
  isw   = getopt(HAM.info,'sweep','isw',  []);
  nsw   = getopt(HAM.info,'sweep','nsw',  []);
  rtol  = getopt(HAM.info,'sweep','rtol', []);

  dloc=getDimQS(HAM.oez(1).op); dloc=dloc(1);

  if isempty(isw) wbdie('missing info isw !?');
  else
     if ~isempty(nsw)
          s{2}=sprintf(' %g/%g',isw,nsw);
     else s{2}=sprintf(' %g',isw); end
     s={[s{:}]};
  end
  if ~isempty(Nkeep), s{end+1}=sprintf('; Nkeep=%g',Nkeep); end
  if ~isempty(rtol ), s{end+1}=sprintf('; rtol=%g', rtol ); end

  q=load_dmrg_info(HAM,'NPsi');
  if q, s{end+1}=sprintf(' keeping %g global multiplets',q);
  else  s{end+1}=''; end

  banner('%s (L=%g%s)%s',s{1},L,['',s{2:end-1}],s{end});

  bupd=any(cellfun(@(x) isequal(x,'-b'), o2site));
  if bupd, o={'-b'}; else o={}; end

  wblog_iter('header',o{:});
  X1=load_dmrg_data(HAM,k1);
  dir=+1;

  for k=k1:L-1
     HAM.info.kc=k;
     X2=load_dmrg_data(HAM,k+1);

     if getqdir(X1.AK,2)<0
        wblog(' * ','search for current site (skip %02g>>)',k);
        if k<L-1, X1=X2; end
        skipped=skipped+1; continue;
     elseif skipped
        X3=load_dmrg_data(HAM,k+2);
        if     X1.info.isw>X3.info.isw, q=+1;
        elseif X1.info.isw<X3.info.isw, q=-1; else q=0; end

        if q>0 || (q==0 && X1.info.time>X3.info.time)
           k1=k; skipped=0; clear X3
        elseif q<0 || (q==0 && X1.info.time<X3.info.time)
           k1=1; k2=k+1; clear X3
           break;
        else
           wbdie('failed to determine direction to continue from kc=%g',k+1);
        end
     end

     [X1,X2,r2,I]=update_psi_2site(HAM,X1,X2,k,k+1,dir,o2site{:});
     RR(k,1)=X2.info.Rho(end);

     if dloc>1 && k<L,
     RR(k+1,3)=X2.info.rho; end

     if ~isempty(X1.HK), HK=getBlockHK(X1);
        if HK
           [eh,Ih]=eigQS(HK); eh=min(eh(:,1));
           Ek=QSpace(Ih.EK)-eh; Ek.info.e0=eh;
           EE(k,1)=Ek;
        end
     end

     save_dmrg_data(HAM,X1,k);
     save_dmrg_data(HAM,X2,k+1);

     if k<L-1, X1=X2; end

     Il(k,1)=wblog_iter(k,dir,L,r2,I,isw);
  end
  wblog_iter('sep');

  if k==L-1
  updateHK(HAM,L,'>>',X2); end

  q=X1; X1=X2; X2=q;

  dir=-1;
  for k=k2:-1:2, if k<L
     HAM.info.kc=k;
     X2=load_dmrg_data(HAM,k-1); end

     [X2,X1,r2,I]=update_psi_2site(HAM,X2,X1,k-1,k,dir,o2site{:});
     RR(k-1,2)=X2.info.Rho(end);

     if dloc>1,
     RR(k-1,4)=X2.info.rho; end

     if ~isempty(X1.HK), HK=getBlockHK(X1);
        if HK
           [eh,Ih]=eigQS(HK); eh=min(eh(:,1));
           Ek=QSpace(Ih.EK)-eh; Ek.info.e0=eh;
           EE(k-1,2)=Ek;
        end
     end

     save_dmrg_data(HAM,X1,k);
     save_dmrg_data(HAM,X2,k-1);

     if k>2, X1=X2; end

     Il(k-1,2)=wblog_iter(k,dir,L,r2,I,isw);
  end

  E0=getdatafield(Il,'E0','-q');
  r2=getdatafield(Il,'r2','-q');

  [~,~]=updateHK(HAM,1,'<<',X2);

  global gHSS
  HAM.info.HSS{isw}=gHSS;

  addto_dmrg_info(HAM,HAM,Il,E0,r2);

  r2=max(r2(find(~isnan(r2))));

  E0=E0(find(~isnan(E0))); E0=[mean(E0), std(E0)];

  wblog_iter('sep');

  tst_matching_bond_spaces(HAM);

end

% -------------------------------------------------------------------- %
function [E,e0]=getEnergy(HAM,k,dir)

   X=load_dmrg_data(HAM,k); odir=itags2odir(X.AK);
   if odir, X.AK.info.itags
      wbdie('got other than current site !? [k=%g: %g]',k,odir);
   end

   L=numel(HAM.mpo);

   if dir>0, k1=k; else k1=k-1; end
   q=getitags(X.AK);

   if dir>0, k1=k; q=regexprep(q{2},'(\d+)\**$','$1a');
      E=QSpace(permuteQS(getIdentityQS(X.AK,1,X.AK,3,q),[1 3 2]));
      Psi=QSpace(contractQS(X.AK,'13',E,'13*'));

      X.AK=E;

      [X1,e0]=updateHK(HAM,k,'>>',X);
      X2=load_dmrg_data(HAM,k+1);

   else k1=k-1; q=regexprep(q{1},'(\d+)\**$','$1b')
      E=QSpace(permuteQS(getIdentityQS(X.AK,2,X.AK,3,q),[3 1 2]));
      Psi=QSpace(contractQS(E,'23*',X.AK,'23'));

      X.AK=E;

      X1=load_dmrg_data(HAM,k-1);
      [X2,e0]=updateHK(HAM,k,'<<',X);
   end

   q=norm(Psi); if abs(q-1)>1E-12
     s=sprintf('input state Psi not normalized !? (%.3g @ %.3g)',q,abs(q-1));
     if abs(q-1)>1E-8, wbdie('%s',s); 
     else wblog('WRN','%s',s); end
   end

   [J,hconj]=get_bondH(HAM,k1);
   [HPsi,E]=get_HPsi(HAM,Psi,X1,X2,J,hconj);

end

% -------------------------------------------------------------------- %
function [E,e0,r2]=getEnergies(HAM,kc)

  L=numel(HAM.mpo); E=nan(L,2); e0=nan(L,2); r2=nan(L,2);
  o={1E-14,1E4};

  for k=kc:L-1
      [E(k,1),e0(k,1)]=getEnergy(HAM,k,+1);
      r2(k,1)=update2Site(HAM,k,k+1,'>>',o{:});
      wblog('  *','(1/3) %3g/%g: e0=%.8g',k,L,E(k,1)/L);
  end

  for k=L:-1:2
      [E(k,2),e0(k,2)]=getEnergy(HAM,k,-1);
      r2(k,2)=update2Site(HAM,k-1,k,'<<',o{:});
      wblog('  *','(2/3) %3g/%g: e0=%.8g',k,L,E(k,2)/L);
  end

  for k=1:kc
      [E(k,1),e0(k,1)]=getEnergy(HAM,k,+1); if k<kc
      r2(k,1)=update2Site(HAM,k,k+1,'>>',o{:}); end
      wblog('  *','(3/3) %3g/%g: e0=%.8g',k,E(k,1)/L);
  end

end

% -------------------------------------------------------------------- %

