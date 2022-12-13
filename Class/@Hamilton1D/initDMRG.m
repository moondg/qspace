function [HAM,Iout]=initDMRG(HAM,varargin)
% function [HAM,Iout]=initDMRG(HAM [,opts])
%
%    initialize A by iDMRG prescription using given
%    Hamiltonian. This initializes both, AK as well as HK.
%    In the end, the orthogonality center is moved from 
%    the center to the left end (kc=1).
%
% Options
%
%    -q            quiet mode
%    -[vV]         verbose mode
%     --complex    enforce complex wave function
%
%    'Nkeep',...   bond dimension to keep (16)
%    'rtol',...    truncation threshold based on density matrix (1E-8)
%    'NPsi',...    number of global states to keep (0)
%    'Qtot',...    global symmetry sector to target
%    'ndav',...    number of Krylov iterations for Davidson (4)
%    'npass',...   number of Davidson basses (4)
%    'nsweep',..   number of block sweeps
%
% Wb,Jul12,20

% wsys='HBL_BaIrO_w1'; J3=[1,.2,.3]; L=16; qloc=1; tflag=1; tst_Hamilton1D

  vflag=0;

  getopt('INIT',varargin);
     Nkeep=getopt('Nkeep',16);
     NPsi=getopt('NPsi', 0);
     Qtot=getopt('Qtot',[]);
     M   =getopt('M',1);

     rtol=getopt('rtol',1E-16);
     ndav=getopt('ndav', 4);
     npass=getopt('npass',4);

     nsw =getopt('nsweep',4);

     if getopt('-v'), vflag=1;
     elseif getopt('-V'), vflag=2; else vflag=0; end

     zflag=getopt('--complex');

  getopt('check_error');

  if size(HAM.ops,2)~=1
     wberr('invalid usage (got %g different sites)',size(HAM.ops,2)); end
  if M<1 || M~=round(M) || numel(M)~=1
     wberr('invalid block size M=%g',M); end

  N=numel(HAM.mpo); M2=2*M;
  if mod(N,M2)
     wberr('invalid length (L=%g @ 2M=%g)',N,M2); end
  L=N/M; L2=L/2;

  o2site=setopts('-',ndav,npass);

  consistency_check(HAM);
  addto_dmrg_info(HAM,HAM);

  E=get_ops_E(HAM,1,'-t'); nq=size(E.Q{1},2);

  lsep={'---',repmat('-',1,50)}; printf('\n'); wblog(lsep{:});
  if NPsi>1
       s=sprintf('%g low-energy states',NPsi);
  else s='ground state'; end
  wblog('<i>','DMRG initalization of %s',s);

  if ~isempty(Qtot)
     if isequal(Qtot,0), Qtot=zeros(1,size(E.Q{1},2));
     elseif ~isnumeric(Qtot) || ~isequal(size(Qtot),[1, nq]), Qtot
        wberr('got invalid Qtot (%g/%g)',numel(Qtot),nq);
     end
     s=[ 'Qtot=[' vec2str(Qtot,'-f') ']' ];
  else s='Qtot=[]';
  end

  if length(s)>20
       wblog(' * ','%s',s); s={'',''};
  else s={[s ', '],''}; end
  if zflag, s{2}='complex'; end

  wblog(' * ','%sNkeep=%g, rtol=%g %s',s{1},Nkeep,rtol,s{2});
  wblog(lsep{:});

  if norm(Qtot) || NPsi>1, Qtot, NPsi, wberr(...
    'non-trivial Qtot or NPsi>1 does not permit wave function prediction!');
  end

  HAM.info.sweep.Nkeep=Nkeep;
  HAM.info.sweep.rtol=rtol;
  HAM.info.sweep.isw=0;

  D1=0;

  for ix=1:L2, l0=(ix-1)*M-1; r0=(L-(ix-1))*M+2;

     HAM.info.sweep.init_DMRG=[l0+(M+2), r0-(M+2)];

     for il=2:M+1, ir=M2-il+3;
        l=l0+il; r=r0-il;

        Xc(il)=load_dmrg_data(HAM,l);
        Xc(ir)=load_dmrg_data(HAM,r);
        if ix>1 || il>2
           Al=getIdentity(Xc(il-1).AK,2,E,[1 3 2]);
           Ar=getIdentity(Xc(ir+1).AK,1,E,[3 1 2]);
        else
           Al=getAtensorLoc(E,[],'L');
           Ar=getAtensorLoc(E,[],'R');
        end

        q=getDimQS(Al);
        if q(1,2)>Nkeep
           if ~D1, D1=q(1,2);
              wblog(' * ',['iDMRG truncating within %s block ' ... 
                '(M=%g, D=%g->%g)'],int2str2(ix,'-th'),M,D1,Nkeep);
              fprintf(1,['\n%8s  ix/L2  isw  site    Nkeep     ' ...
                'E0        e0      drho2     exp(SE)\n']);
           else D1=max(D1,q(1,2)); end

           Q=get_center_randn(Al,Ar,zflag);
           [U,S,Vd]=svdQS(Q,2,'Nkeep',Nkeep);
           Al=contract(Al,2,U,1,[1 3 2]);
           Ar=contract(Vd,2,Ar,1);
        end

        Xc(il).AK=setitags(Al,'-A',l);
        Xc(ir).AK=setitags(Ar,'-A',r);

        Xc(il)=updateHK(HAM,l,'>>',Xc(il));
        Xc(ir)=updateHK(HAM,r,'<<',Xc(ir));

        save_dmrg_data(HAM,Xc(il),l);
        save_dmrg_data(HAM,Xc(ir),r);
     end

     if D1 || ix==L2
        Q=get_center_randn(Xc(il).AK,Xc(ir).AK,zflag);
        Ar=contract(Q,2,Xc(ir).AK,1);

        q=normQS(Ar); if q<1E-8
          wbdie('got |Xc(%g).AK| = %.2g (l=%g, M=%g)',ir,q,l,M); end
        Xc(ir).AK=Ar;

        nsw2_1=2*nsw+1;
        cv=ones(1,M2+2);
        if ix>1, i2=[2 M2]; else i2=[3 M2-1]; end

        for isw=1:nsw2_1, sdir=(-1)^(isw-1);
           if isw==1, iL=M+1:i2(2);
           elseif isw==nsw2_1, iL=il-1:M+1;
           elseif sdir>0, iL=il-1:i2(2);
           else iL=il+1:-1:i2(1); end

           for il=iL, l=l0+il; o={};
              if il==M+1, r=r0-il; o={'-i'};
              elseif il<=M
                 r=l+1;
              else
                 r=r0-(M2+2-il);
                 l=r-1;
              end
              [Xc(il),Xc(il+1),r2,I2]=update_psi_2site(HAM,...
               Xc(il),Xc(il+1),l,r,sdir,o{:},o2site{:});

               save_dmrg_data(HAM,Xc(il  ),l);
               save_dmrg_data(HAM,Xc(il+1),r);
              cv(il)=I2.converged;
           end

           s={'',''};
           if sdir>0
                s{1}=sprintf('|%02g>',l);
           else s{1}=sprintf('<%02g|',l); end

           d=getDimQS(Xc(il).AK); if size(d,1)>1
                s{2}=sprintf('%g/%g',d(:,2));
           else s{2}=sprintf('%g',d(:,2)); end
           E0=I2.E0(end); E0(2)=E0/(ix*M2);

           if mod(isw,2), fprintf(1,...
             '%s %3g/%-3g %g/%g %5s %8s %10.3f %7.4g  %7.2e  %6.4g\n',...
             datestr(now,'HH:MM:SS'), ix,L2, (isw-1)/2, nsw,...
             s{:},E0,I2.svd2tr,exp(SEntropy(I2.svd)));
           end

           if M==1, break; end
        end
     end

     if ix<L2
        Xc(1)=Xc(il);
        if D1
          [~,Xr,Iu]=update2Site(HAM,Xc(il),Xc(ir),'<<',rtol,Nkeep,'-i');
        else Xr=Xc(ir); end

        q=itags2int(Xr.AK);
        if q(3)~=r, wbdie('unexpected AK data'); end
        Xr.AK=setitags(Xr.AK,'-A',r);
        Xr=updateHK(HAM,r,'<<',Xr);

        Xc(M2+2)=Xr;

        save_dmrg_data(HAM,Xr,r);
     end
  end

  Xl=Xc(il); Xr=Xc(ir);
  for l=r-1:-1:1
     [Xl,Xr,Iu]=update2Site(HAM,Xl,Xr,'<<',rtol,Nkeep);
     save_dmrg_data(HAM,Xr,l+1);
     if l>1
          Xr=Xl; Xl=load_dmrg_data(HAM,l-1);
     else save_dmrg_data(HAM,Xl,l);
     end
  end

  HAM.info.sweep.init_DMRG=1;

  addto_dmrg_info(HAM,NPsi);

  Iout=add2struct('-',Nkeep,rtol,Nkeep,rtol,NPsi,Qtot,zflag,ndav,npass);

end

% -------------------------------------------------------------------- %
function Xk=fix_HK(Xk,k,dir)

   if ~isnumeric(dir), dir=check_dir(dir); end

   da=getDimQS(Xk.AK); da=da(1,:);
   dh=getDimQS(Xk.HK); if ~isempty(dh), dh=dh(1,:);
   else
      wblog('WRN','got empty H(%g)',k);
      dh=0;
   end

   if dir>0
      if any(dh~=da(2))
         if itags2odir(Xk.AK)~=2
            wberr('invalid A-tensor'); 
         end
         wblog(' * ','adding zero-blocks to HK(%g): D=%g->%g',k,dh(1),da(2));
         E=QSpace(contractQS(Xk.AK,'13*',Xk.AK,'13'));
         Xk.HK=Xk.HK+1E-24*E;
      end
   else
      if any(dh~=da(1))
         if itags2odir(Xk.AK)~=1
            wberr('invalid A-tensor'); 
         end
         wblog(' * ','adding zero-blocks to HK(%g): D=%g->%g',k,dh(1),da(1));
         E=QSpace(contractQS(Xk.AK,'23*',Xk.AK,'23'));
         Xk.HK=Xk.HK+1E-24*E;
      end
   end

end

% -------------------------------------------------------------------- %

function X=get_center_randn(Al,Ar,zflag)

   X=getIdentityQS(Al,2,'-0');
   [i1,i2]=matchIndex(X.Q{2},Ar.Q{1});
   for i=1:numel(i1)
      s=[ size(X.data{i1(i)},2), size(Ar.data{i2(i)},1) ];
      X.data{i1(i)}=zeros(s);
   end

   if zflag, o={'-z'}; else o={}; end
   X=randomize(QSpace(X), 1., o{:});

   X=setitags(X,2,Ar,1);
end

% -------------------------------------------------------------------- %

