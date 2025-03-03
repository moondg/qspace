function [Htot,Iout]=runED(HAM,varargin)
% function [Htot,Iout]=runED(HAM [,opts])
%
%    Perform exact diagionalization (ED) on given Hamiltonian.
%
%    NB! ED here is intended for the full diagonalization.
%    if only a particular set of lowest energy states shall be
%    targeted since system is too large for full diaognalization,
%    use NPsi>1 with regular DMRG run.
% 
% Options
% 
%    -k       stop with keyboard at end of routine
%    -v       verbose mode (report progress)
% 
% Wb,Dec28,21 ; Wb,Mar03,22

% for last tentative version to also perform iterative Krylov
% see Archive/runED_220303.m // Wb,Mar03,22

% wsys='Kitaev'; L=2; W=4; J2=0.5; useAC=1; tflag=2; tst_Hamilton1D
% [E,Iout]=runED(HAM,'-k','npass',12,'NPsi',32);

   getopt('INIT',varargin);
      if getopt('-v'), vflag=2;
      elseif getopt('-q'), vflag=0; else vflag=1; end
      pflag=getopt('-p');

      Dmax=getopt('Dmax',1E5);
      fflag=getopt('-F');

      kflag=getopt('-k');

   getopt('check_error');
   L=length(HAM);

   if size(HAM.oez,2)>1, wbdie(...
     'invalid usage (setting with different sites not yet implemented)'); end
   if isstruct(HAM.mpo(1)), wbdie('got pseudo MPO (full MPO required)'); end
   hconj=HAM.info.mpo.use_hconj; % wheth H_full = H + H'

   E=HAM.oez(1).op;
   dloc=getDimQS(E); dloc=max(dloc,[],2); D=dloc.^L;
   if D(1)>Dmax || D(end)>Dmax && ~fflag, wbdie(...
     'Hilbert space too large [ (d=%d)^%d = %.3gG ]',dloc(end),L,D(end)/2^30);
   end

   L2=ceil(L/2);

   if vflag
   wblog(' * ','building basis and matrix elements (L=%g)',L); end

   for l=1:L2,
      if l==1, Al=getvac(E); end
      Al=getIdentity(Al,2,E,[1 3 2]);

      if vflag, q=getDimQS(Al); D=q(1,2);
         fprintf(1,'\r  %2d/%d D=%4g \r',l,L2,D); end

      Xl=load_dmrg_data(HAM,l);
         Xl.AK=setitags(Al,'-A',l); if l==L2
         Xl.AK=itagrep(Xl.AK,2,'^K','L'); end
      Xl=updateHK(HAM,l,'>>',Xl);
      save_dmrg_data(HAM,Xl,l);

      r=L-l+1; if r<=L2, continue; end

      Xr=load_dmrg_data(HAM,r);
         Xr.AK=setitags(permute(Al,[2 1 3]),'-A',r); if r==L2+1
         Xr.AK=itagrep(Xr.AK,1,'^K','R'); end
      Xr=updateHK(HAM,r,'<<',Xr);
      save_dmrg_data(HAM,Xr,r);
   end

   if vflag
   wblog(' * ','building full system unitary (D=%d)',dloc(end)^L); end

   ALR=getIdentity(Xl.AK,2,Xr.AK,1,[1 3 2]);
   [Qtot,DD,dc]=getQDimQS(ALR,2);

   Htot=getIdentity(ALR,2); Htot.data=repmat({0},size(Htot.data));
   if ~isequal(Htot.Q{1},Qtot), wbdie('Htot not sorted in QIDX !?'); end
   if vflag, qfmt=['(' getqfmt(Htot) ')']; end

   nQ=size(Qtot,1);
   for k=1:nQ
      ULR=getsub(ALR, matchIndex(ALR.Q{2},Qtot(k,:)));
      if vflag, q=getDimQS(ULR);
         if size(q,1)>1
              s=sprintf('D* = %5d / %5d (%5d)',DD(k),q([1 end],2));
         else s=sprintf('D= %5d / %5d',DD(k),q(2)); end
         fprintf(1,['   %2d/%d q=' qfmt '   %s\n'],k,nQ,Qtot(k,:),s);
      end
      H=contract(ULR,'!2*',{Xl.HK,{ULR,Xr.HK}});

      if isempty(H), H.data{1}=zeros(1,DD(k)); continue
      elseif numel(H.data)~=1, H, wbdie('unexpected H'); end

      H=H.data{1};
      if hconj
         Htot.data{k}=eig(H+H')';
      else
         e=norm(H-H','fro')/norm(H,'fro');
         if e>1E-8, wbdie('H not hermitian @ e=%.3g',e); 
         elseif e>1E-12, wblog('WRN','H not hermitian @ e=%.3g',e); end
         Htot.data{k}=eig(H+H')'/2;
      end
   end

   ee=Htot.data; ee=sort([ee{:}])';

   deps=max(1E-12,1E-8*(ee(end)-ee(1))); D=numel(ee);
   i=[0; find(diff(ee)>deps); D];
   di=diff(i); gmax=max(di); e=ones(size(di));
   gg=sparse(e,di,e,1,gmax);

   if gmax>8
      wblog('NB!','got degeneracy %g (using eps=%.3g)',gmax,deps);
   elseif gmax>1, wblog(' * ','got degeneracy <= %g (eps=%.2g)',gmax,deps); 
   else
      de=min(diff(ee));
      s={ sprintf('non-degenerate spectrum @ %.3g',de)
          sprintf('eps=%.1g',deps) };
      if ~isempty(dc), m=max(prod(dc,2)); else m=1; end
      if m>1, wblog(' * ','%s (m<=%g; %s)',s{1},m,s{2});
      else    wblog(' * ','%s (%s)',s{:}); end
   end

   if nargout>1, Iout=add2struct('-',ee,gg,gmax,deps,di); end

   if ~pflag
      if kflag, keyboard, end
      return
   end

   p=HAM.info.param;
   if isfield(p,'istr'), s=untex(p.istr); else s=''; end

ah=smaxis(2,2,'tag',mfilename); addt2fig Wb
header('%M :: %s',s);
header('SW',param2str(p,'-x','system'),{'FontSize',10});

   setax(ah(1,1)); j=100;
   if D>2*j
      plot(1:j,ee(1:j),'.-'); grid on
      setax(ah(1,2));
   end
   plot(ee,'.-'); xtight(1.1); ytight(1.1); grid on
   title2('full spectrum (D=%g)',D);

   if gmax>1
   setax(ah(2,1)); 
      j=find(di==max(di),1);
      ig=i(j)+1:i(j+1);
      plot(ig,ee(ig),'o-');
      title2('largest degeneracy @ %g',gmax);
   end

   setax(ah(2,2));
   plot(full(gg))
   label('degeneracy','# occurances','histogram on degeneracies');

   if nargout>1, Iout.ah=ah; end
   if kflag, keyboard, end
end

% -------------------------------------------------------------------- %

