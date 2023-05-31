function [X1,X2,r2,Iout]=update_psi_2site(HAM,X1,X2,k1,k2,kdir,varargin)
% function [X1,X2,r2,Iout]=update_psi_2site(HAM,X1,X2,k1,k2,kdir [,opts])
%
%    DMRG 2-site update step by switching to an enlarged bond picture,
%    If option -b if set, then a plain, ie., non-enlarged, bond update
%    is perfomed, instead.
%
% Wb,Apr10,14

% see also Hamilton1D/runED.m
% which adapted the Davidson algorithm below // Wb,Dec30,21

% -------------------------------------------------------------------- %
% also allow global state index Psi with current site // Wb,Jul28,16
% -------------------------------------------------------------------- %
% outsourced from sweepPsi.m // Wb,Jan20,16:
% for old version, see Hamilton1D/Archive/sweepPsi_160119.m
% see also $MLAB/david.m and
% david2SiteQS.cc for former C-version for abelian symmetries
% -------------------------------------------------------------------- %

  rtol  = getopt(HAM.info,'sweep','rtol', 1E-15);
  Nkeep = getopt(HAM.info,'sweep','Nkeep',   32);
  isw   = getopt(HAM.info,'sweep','isw',     []);
  vflag=1;

  getopt('init',varargin);
     if getopt('-q'), vflag=0;
     elseif getopt('-v'), vflag=2; end

     npass = getopt('npass', 1);
     ndav  = getopt('ndav',  4);
     rtol  = getopt('rtol', rtol);
     iflag = getopt('-i');
     bupd  = getopt('-b');

  getopt('check_error');

  tcpu=cputime;

  if ~iflag
     q=[ itags2odir(X1.AK), itags2odir(X2.AK) ];
     if isempty(matchIndex(q,[2 0; 0 1])), save2tmp -5
        X1.info.itags, X2.info.itags, q, wbdie(...
       'ERR got invalid orthonormalization [ %g, %g ]', q(1),q(2));
     end
     if k2~=k1+1, wbdie(['invalid usage ' ... 
       '(expecting two consecutive sites [%g,%g])'],k1,k2);
     end
  end

  L=numel(HAM.mpo);

  if kdir>0
       d=getDimQS(X1.AK); d=d(end,:); d([1 2])=d([2 1]);
  else d=getDimQS(X2.AK); d=d(end,:); end

  if rtol<0
     warning('Wb:WRN','got rtol=%g !? (set to 0)',rtol);
     rtol=0;
  end

  i=[ d(1)>=prod(d(2:end))
      d(2)>=prod(d([1,3:end]))
  ];

  if 0 & any(i)
   % WRN this can lead to trouble if Nkeep is increased
   % at once by more than a factor of d_local!
   % => this would never update below!
   % -------------------------------------------------------------- %
   % if i(2) // don't (ortho2site() assumes *current* site)
   %  % make sure got minimum bond dimension
   %    [X1.AK,X2.AK,r2,Iout]=ortho2site(X1.AK,X2.AK,-kdir,rtol,Nkeep);
   % end

     [X1.AK,X2.AK,r2,Iout]=ortho2site(X1.AK,X2.AK, kdir,rtol,Nkeep);
     if kdir>0
          X1=updateHK(HAM,k1,'>>',X1); Iout.ex=X1.info.ex;  
     else X2=updateHK(HAM,k2,'<<',X2); Iout.ex=X2.info.ex;  end

     return
  end

  tpat='[a-zA-Z]+(\d+)';

  d1=itags2odir(X1.AK);
  d2=itags2odir(X2.AK);
  if d1 && d2, wbdie('got current site other than A1 and A2 !?'); end

  Ek1=get_local_id(HAM,k1,X1.AK,3);
  Ek2=get_local_id(HAM,k2,X2.AK,3);

  if ndav>0 && ~bupd
     t=regexprep(getitags(X1.AK,2),tpat,'E$1l');
     E1=QSpace(permuteQS(getIdentityQS(X1.AK,1,Ek1,t),[1 3 2]));
     A1=QSpace(contractQS(E1,'*',X1.AK));

     X1.AK=E1;
     X1=updateHK(HAM,k1,'>>',X1);

     t=regexprep(getitags(X2.AK,1),tpat,'E$1r');
     E2=QSpace(permuteQS(getIdentityQS(X2.AK,2,Ek2,t),[3 1 2]));
     A2=QSpace(contractQS(E2,'*',X2.AK));

     X2.AK=E2;
     X2=updateHK(HAM,k2,'<<',X2);

     r1=numel(A1.Q); r2=numel(A2.Q);
     NPsi=(r1>2 || r2>2);
     if r1<2 || r2<2 || NPsi && r1+r2~=5
        wbdie('invalid usage (got [%g,%g]-rank A-tensors)',r1,r2); 
     end
     p=[]; if r1==3, p=[1 3 2]; end

     Psi=contract(A1,A2,p);
  else
     if bupd, o={'-b'}; else o={}; end

     [X1.AK,X2.AK,r2,I2]=ortho2site(X1.AK,X2.AK,kdir,rtol,Nkeep,o{:});
     if kdir>0 || bupd, X1=updateHK(HAM,k1,'>>',X1); end
     if kdir<0 || bupd, X2=updateHK(HAM,k2,'<<',X2); end

     if kdir>0
          Iout.ex=X1.info.ex;
     else Iout.ex=X2.info.ex; end

     if ndav<=0, return; end

     Psi=I2.S;
  end

  dpsi=getDimQS(Psi); rpsi=size(dpsi,2);
  if rpsi>2
       oc='!3*'; NPsi=dpsi(1,end); NPsi2=dpsi(end);
  else oc='*';   NPsi=0; NPsi2=1; end
  NPsi1=max(1,NPsi);

  if NPsi
       G0=QSpace(getIdentityQS(Psi,3));
  else G0=[];
  end

  R=contract(Psi,oc,Psi);
  if rpsi>2
       e=normQS(getIdentityQS(R)-R);
  else E=blkdiag(R); e=norm(E-eye(size(E)));
  end

  reps=1E-12*prod(dpsi(1,:))^(1/rpsi);

  stol=sqrt(rtol);

  if e>reps
     s=sprintf('Psi(%g:%g) not normalized @ %.3g / %.3g',k1,k2,e,reps);

     if e>1E-6
        nrm2=normQS(Psi)^2 / NPsi2;
        s={ s, sprintf('|Psi|^2 = %.3g @ NPsi=%d',nrm2,NPsi2) };
        if e>0.01
           banner(['ERR ',s{1}]);
           m=sprintf('./tmp_%s_%g_%g_%s_checkNorm.mat',...
             mfilename,k1,k2,iff(kdir>0,'lr','rl'));
           save2(m,'-f');
           wbdie('%s !?',s{2});
        elseif e>1E-6, wblog('WRN %s\n%s',s{:});
        end

     else wblog('WRN',s); end
  end

  E0=[]; H=[]; Ig=[]; converged=0; Psi0=Psi;

  pseudo_mpo=(~usingFullMPO(HAM));
  if pseudo_mpo, [J,hconj]=get_bondH(HAM,k1);
       oH={J,hconj};
  else oH={'--full-mpo'};
     q=getfield2(HAM.info,'sweep','init_DMRG',{0});
     if numel(q)==2 && all(q>1)
        oH{1}=[oH{1},'-INIT'];
     end
  end

  if ndav>=prod(dpsi(1,:))
     d1=getDimQS(E1);
     d2=getDimQS(E2);
     if vflag>1, wblog('TST',...
        'david(%g,%g) with ndav=%g @ d=[%s; %s] {%g,%g}',k1,k2,ndav,...
         vec2str(dpsi(1,:),'-f'), vec2str(dpsi(2,:),'-f'),d1(1,2),d2(1,1));
     end
     q=max([d1(1,2); d2(1,1)]); if ndav>q, ndav=q; end
  end

  if ~iflag && npass>10, wblog('WRN','got npass=%g !?',npass); end
  if ndav<1 || ndav>10, wblog('WRN','got ndav=%g !?',ndav ); end

  if npass<1
     [HPsi,E0]=get_HPsi(HAM,Psi,X1,X2,oH{:});
     E0=[E0,E0];
  end

  for ipass=1:npass, Ak=Psi;
     if NPsi
        setitags(Ak,'dav0',3);
     end

     for idav=1:ndav
      % need to calculate H|psi> twice (see HAM1 and HAM2 below)
      % (1) given input or improved ground state Psi from (2),
      %     find new search direction, i.e. n -> n+1.
      % (2) given new search direction |psi>:=|n>:=|A(idav)>,
      %     calculate <n'|H|n> to grow H by one row and column
      %  => Psi := improved new ground state from enlarged H
      %
      % This therefore differs from Krylov [even in the absence
      % of the conditioning step below since Krylov would include
      % both H|psi> iterations to grow H; nevertheless the result
      % should be equivalent when using a Krylov subspace of
      % dimension 2*ndav since after all what is used is H^n|psi>].
      % See also $MLAB/david.m // Wb,Aug10,16

      % ---------------------------------------- %
        HPsi=get_HPsi(HAM,Psi,X1,X2,oH{:});
        E=contract(Psi,oc,HPsi);

        q=blkdiag(E); e=norm(q-q')/max(1,norm(q));
        if e
           if e>1E-12
              wbdie('got non-Hermitian H (e=%.3g) !?',e); end
           if NPsi, E=0.5*(E+E'); else E=real(getscalar(E)); end
        else
           if ~NPsi, E=getscalar(E); end
        end

        if idav==1
           if NPsi, HD=struct(E);
              for i=1:numel(HD.data), HD.data{i}={HD.data{i}}; end
           else
              HD=struct(getvac(Psi)); HD.data={{E}};
           end
        end

      % find new search direction (Davidson algorithm)
      % from currently best known approximation to ground state
      % NB! here without conditioning, i.e.(E-diag(H))\[...] is
      % not included, since in the presence of non-abelian symmetries,
      % "diag(H)" is not simply acessible for the coupling terms.
      % see als Archive/update_psi_2site_160810.m for older version.
      % Wb,Apr12,14 ; Wb,Aug11,16

        if ~NPsi
             Q=HPsi-E*Psi;
        else Q=HPsi-contract(Psi,E);
        end

      % NB! certain symmetry sectors may have (1) converged,
      % or (2) exhausted their available state space, while
      % others still may proceed => need selective truncation
      % of new search space via orthonormalization
      % (see also 2nd orthonormalization below)
        q=norm(E)/sqrt(NPsi2); x=0;
        [Q,ns,nrm]=ortho_Psi3(Q,q*stol,NPsi);

        if isempty(Q)
           if vflag
              if NPsi, s='all states'; else s='state'; end
              fprintf(1,'> %s  %s converged (%d/%d, %d/%d: %.2g/%.2g)  \n',...
              datestr(now,'HH:MM:SS'),s,ipass,npass,idav,ndav,nrm*nrm,q*reps);
           end

           if ~NPsi, Ek=E; else Ek=eigQS(E); end
           H1=Ek; converged=1;
           break
        else
           for p=1:2
              for j=1:idav
                 Q=ortho_PsiX(Q,Ak(j),NPsi);
                 if isempty(Q.data), break; end
              end

              if p>1 && normQS(Q)^2<NPsi2*reps
                 x(p)=0;
              else
                 [Q,ns,nrm]=ortho_Psi3(Q,sqrt(NPsi2)*stol,NPsi);
                 if ns, x(p)=nrm/sqrt(ns); else x(p)=0; end
              end
              if x(p)==0, break; end
           end
        end

        if norm(x)<1E-12
         % already got sufficiently tiny number in p==1 pass
         % e.g. this can happen towards the ends of the chain where
         % the Hamiltonian has degeneracies (such as "zero-blocks")
         % which cannot be further resolved  by simply applying H
         % to given state space // Wb,Apr15,14
           if vflag>1, d=getDimQS(Psi); wblog('WRN',[...
              'david() got insufficient state space\n' ... 
              'having idav=%g/%g @ %.3g; d=[%s; %s]'], idav, ndav, ...
              norm(x),vec2str(d(1,:),'-f'),vec2str(d(end,:),'-f'));
           end
           ndav=idav-1; break;
        elseif ~x(p)
           wblog('KRL','state space fully exhausted (x=%.3g)',norm(x));
           ndav=idav-1; break;
        elseif abs(x(p)-1)>1E-2, wblog('WRN',...
          'got non-orthonormal Krylov space (e=%.3g) !?',abs(x(p)-1));
        end

        if NPsi
           setitags(Q,sprintf('dav%g',idav),3);
        end
        l=idav+1; Ak(l)=Q;

        HPsi=get_HPsi(HAM,Ak(l),X1,X2,oH{:});

        for j=1:l
           h=contract(Ak(j),oc,HPsi);
           HD=setHD_block(HD,j,l,h);
        end

        [H,HD]=fuse_data(HD);

        if NPsi<=1
           if isnumeric(H), Hk=H;
           else
              if numel(H.data)~=1 || ~isequal(size(H.data{1}),[l l])
                 wbdie('unexpected H object !?'); end
              Hk=H.data{1};
           end
           H1=min(diag(Hk));

           [Uk,Ek]=eig(Hk);

           Psi=QSpace;
           for j=1:l
              Psi=Psi+Uk(j,1)*Ak(j);
              if j==1 && NPsi2>1
                 Psi=setitags(Psi,'',3);
              end
           end
           if NPsi2>1
              Psi=setitags(Psi,3,Psi0,3);
           end
        else
           [ee,Ig]=eigQS_dav(H,NPsi1,G0);

           if Ig.NK~=NPsi1
              warning('Wb:WRN','kept %g/%g multiplets !?',Ig.NK,NPsi1); end
           Ek=sort([Ig.EK.data{:}]);
           H1=min(diag(blkdiag(H)));

           UD=split_UDAV(Ig.AK,HD);

           Psi=QSpace;
           for j=1:l, Psi=Psi+contract(Ak(j),UD(j)); end
        end

        if converged || idav>ndav, break; end

     end

     E0(ipass,:)=[Ek(1),H1(1)];
     if converged, break; end

  end

  if NPsi, d=getDimQS(Psi);
     if size(d,2)~=3, Psi, d, wbdie('invalid NPsi !?'); end
     if d(end)~=NPsi2
        wblog('NB!',...
          'low-energy symmetry multiplets changed (%d->%d, %d->%d)',...
           NPsi1,d(1,end),NPsi2,d(end));
        display(R(1)); display(QSpace(getIdentityQS(Psi,3)));

        NPsi1=d(1,end); NPsi2=d(end); NPsi=NPsi1; 
     end
  end

  oe={'stol',stol};
  if Nkeep>0, oe(3:4)={'Nkeep',Nkeep}; end

  [U,B,Iout]=orthoQS(Psi,iff(kdir>0,1,2),oe{:});

  if NPsi && Iout.svd2tr>1E-32
     B=ortho_Psi3(B,stol^2,NPsi1);
  end

  if NPsi>1
     if isempty(Ig)
        [H,HD]=fuse_data(HD);
        [ee,Ig]=eigQS_dav(H,NPsi1,G0);
     end

     Eg=QSpace(Ig.EK);
     Eg.info.itags=regexprep(Eg.info.itags,'K(\**)$','$1');

  else Eg=Ek(1); end

  if converged, converged=idav; end
  add2struct(Iout,Eg,converged);

  if 0, clear R E
     for i=1:2
        if i==1, q=Psi; s='Psi'; else q=QSpace(B); s='B'; end
        R(i)=contract(q,oc,q);
        E{i}=blkdiag(R(i)); e=norm(E{i}-eye(size(E{i})));

        s=sprintf('got orthonormal %-3s @ %.3g',s,e);
        if e>1E-8
             warning('Wb:WRN',[s ' !?']);
        else wblog('TST',[s ' (ok)']); end
     end
  end

  if size(Iout.svd,2)==2
     rr=[ Iout.svd(:,1).^2 / NPsi2, Iout.svd(:,2) ];
     r1=rr(:,1).*rr(:,2); i=find(r1>0);
     se=rr(i,2)'*(-rr(i,1).*log(rr(i,1)));
  else
     rr=Iout.svd.^2 / NPsi2; r1=rr; i=find(r1>0);
     se=-rr(i)'*log(rr(i));
  end

  e=sum(r1); e(2)=1-e; e(3)=abs(e(2));
  if e(3)>1E-8
     if e(3)>=1E-3, save2('./tmp.mat'); end
     s=sprintf('%.8g @ %.3g',e(1),e(2));
     q=getfield2(HAM.info,'sweep','isw',{0});
     if e(3)<1E-3 || q==1 && e(3)<1
          wblog('WRN','wave function normalized @ %s',s);
     else wbdie('%s wave function not normalized (%s)',s);
     end
  end

  Iout.rr=r1; i=find(r1>0);
  Iout.se=se;

% NB! in the adaptive approach above using eigQS(),
% keeping the same number of NPsi multiplets can lead
% to a different number of states as compared to NPsi2.
% NB! checking integer e also ensures that the Davidson
% basis A above is othonormal, indeed.

  q=sum(r1); e=abs(q-round(q));
  [q(2),q(3)]=get_NPsi_bond(B);

  if e>1E-3 || q(2)~=NPsi, warning('Wb:WRN',...
    'got invalid normalization sum(rr)=%g !?',sum(r1)/NPsi2); end
  if q(3)~=NPsi2, wblog('PSI',...
    'having %g multiplets, yet %g->%g states',NPsi1,NPsi2,q(3));
  end

  if nargout>3, Iout=add2struct(Iout,kdir,H,E0); end

  r2=Iout.svd2tr;
  dsw=[ Nkeep, Iout.se, [E0(end,1) std(E0(:))]/L, r2 ];

  if kdir>0
     if NPsi, p=[1 3 4 2]; else p=[]; end
     X1.AK=contract(X1.AK,2,U,1,[1 3 2]);
     X2.AK=contract(B,2,X2.AK,1,p);

     if pseudo_mpo
          X2.OP=QSpace(size(X2.OP));
     else X2.HK=QSpace; end

  elseif kdir<0
     if NPsi, p=[1 3 2 4]; else p=[1 3 2]; end
     X2.AK=contract(U,1,X2.AK,1);
     X1.AK=contract(X1.AK,2,B,1,p);

     if pseudo_mpo
          X1.OP=QSpace(size(X1.OP));
     else X1.HK=QSpace; end
  end

  pat='^[A-Z]+(\d+)[A-Za-z]*';
  X1.AK=itagrep(X1.AK,2,pat,'K$1');
  X2.AK=itagrep(X2.AK,1,pat,'K$1');

  if check_Qmatch(X1.AK,2,X2.AK,1), save2(...
     sprintf('./tmp_%02g_%02g_%s1_Qmatch',k1,k2,iff(kdir>0,'lr','rl')));
  end

  sc={'!1','!2','!3'}; Rfac=0;
  if numel(Psi.Q)>2 && isempty(regexp(Psi.info.itags{end},'\d')) % `Psi'
     q=getDimQS(Psi); q=q(:,end);
     if q(end)>1, sc={'23','13','12'}; Rfac=1/q(end); end
  end

% add reduced density matrix 'rho' for local state space
% which also allows to double check whether local state space
% has dropped out along the simulation! // Wb,Mar26,15
% NB! AK is (significantly) smaller in size (factor 10 for SU(4))
% => rather store AK than expanded bond space % Wb,Mar30,15
  if kdir>0
     DX=getDimQS(Psi); Psi=X2.AK;
     Dk=getDimQS(Psi); Dk=Dk(:,1);

     rho=reduce_rho_q0(Psi,sc{3},Eg);

     Rho=reduce_rho_q0(Psi,sc{1},Eg);
     if Rfac, Rho=Rfac*Rho; end

     I2=X2.info;
        if k2<L, l=2*max([0, isw-1])+1; else l=max([1, isw]); end
        I2.dsw(l,:)=dsw;
        I2.Rho(l,1)=Rho;

        if isempty(I2.Eg), I2.Eg=Eg; end
        I2.Eg(l,1)=Eg;

     I2.itags=X2.AK.info.itags;
     X2.info=add2struct(I2,isw,DX,Dk,rho,rr,r2);
     X2.Psi=Psi;
  else
     DX=getDimQS(Psi); Psi=X1.AK;
     Dk=getDimQS(Psi); Dk=Dk(:,2);

     rho=reduce_rho_q0(Psi,sc{3},Eg);
     Rho=reduce_rho_q0(Psi,sc{2},Eg);
     if Rfac, Rho=Rfac*Rho; end

     I1=X1.info;
        if k1>1, l=2*max([0, isw-1])+2; else l=max([1, isw])+1; end
        I1.dsw(l,:)=dsw;
        I1.Rho(l,1)=Rho;

        if isempty(I1.Eg), I1.Eg=Eg; end
        I1.Eg(l,1)=Eg;

     I1.itags=X1.AK.info.itags;
     X1.info=add2struct(I1,isw,DX,Dk,rho,rr,r2);
     X1.Psi=Psi;
  end

  Iout=add2struct(Iout,DX,Dk);

 tcpu(2)=cputime;

 if kdir>0
    [X1,Iout.e0]=updateHK(HAM,k1,kdir,X1,X2); tcpu(3)=cputime;
    X1.info.time=now;
    X1.info.cpu=diff(tcpu); Iout.ex=X1.info.ex;
 else
    [X2,Iout.e0]=updateHK(HAM,k2,kdir,X2,X1); tcpu(3)=cputime;
    X2.info.time=now;
    X2.info.cpu=diff(tcpu); Iout.ex=X2.info.ex;
 end

 if check_Qmatch(X1.AK,2,X2.AK,1), save2(...
    sprintf('./tmp_%02g_%02g_%s2_Qmatch',k1,k2,iff(kdir>0,'lr','rl')));
 end

end

% -------------------------------------------------------------------- %
function e=check_Qmatch(A1,k1,A2,k2)

  [i1,i2,Im]=matchIndex(A1.Q{k1},A2.Q{k2}); e=0;

  if ~isempty(Im.ix1) || ~isempty(Im.ix2)
     w=[0 0]; a1=[]; a2=[]; e=1;

     if ~isempty(Im.ix1), a1=getsub(A1,Im.ix1); w(1)=norm(a1)^2; end
     if ~isempty(Im.ix2), a2=getsub(A2,Im.ix2); w(2)=norm(a2)^2; end

     wblog('WRN',['got dropping out symmetry sector\n ' ... 
       'A#%g: {%s} @ %.8g\n B#%g: {%s} @ %.8g'],...
        k1, vec2str(Im.ix1), w(1), k2, vec2str(Im.ix2), w(2));
     disp(a1), disp(a2)

  end

end

% -------------------------------------------------------------------- %

