function [H0,Iout]=initNRG(HAM,varargin)
% function [H0,Iout]=initNRG(HAM [,opts])
%
%    Initialize A-tensors in DMRG matrix product state by NRG like
%    prescription of iterative diagonalization using given Hamiltonian.
%    This initializes both, AK as well as HK. In order to generate
%    an RL orthonormalized state with kc=1, the procedure starts
%    at the right end of the system (k=L) and stops at k=1. At k=1
%    then the DMRG states to be targeted globally are selected.
%
% Options
%
%   '-v'        more verbose mode
%   'Nkeep',..  number of states/multiplets to keep during NRG iterations (16)
%   'rtol',...  tolerance used with ortho2site()
%   '--maxS'    use fully symmetrized superposition within ground state space
%   '--cplx'    make data complex (even if Hamiltonian is real)
%
% Final targeting of states at iteration / site k=(1 <- L):
%
%   'Qtot',..   (ground state) symmetry sector to choose;
%               if empty, the symmetry sector(s) of the NRG
%               low-energy states is/are taken (default:
%               [] if NPsi or dQtotN is set, all-zero otherwise)
%
%   'NPsi',...  global number of states/multiplets to target
%
%      Here NPsi>0 enforces global Psi index, including NPsi=1.
%      Default: NPsi=0, in which case if also Qtot is not set
%      this tries to initialize the DMRG starting state in
%      the NRG ground state in the global scalar symmetry sector
%      also representing the vacuum state, namely with Qtot all-zeros.
%      In this case no additional index for the global state
%      is required. For NPsi>1, a global index is required
%      in any case as this simultaneously targets multiple states.
%      Its label / itag by default is set to `Psi'. This label
%      is also expected later when performing DMRG sweeps for the
%      case that multiple states are targeted.
%
%   'dQtotN',.. is alternative to NPsi above, specifying [d(Qtot),nPsi]
%
%      dQtotN may contain multiple rows and thus permits to explicitly
%      set the number of multiplets to target in the global symmetry
%      sectors specified. The additional trailing column represents
%      nPsi=dQtotN(:,end) such that NPsi = sum(nPsi); the specified
%      symmetry sectora are relative to Qtot, hence the naming `d(Qtot)'.
%      Same as with NPsi above, if Qtot is not specified, the symmetry
%      sector of the `NRG ground state' is taken for Qtot.
%
%      Note that for NPsi>1 or dQtotN the distribution ov states or
%      multiplets over the symmetry sectors may change via truncation
%      in the Davidson algorithm in the DMRG sweeps. However, the DMRG
%      ensures that in each symmetry sector chosen during initialization
%      at least one state or multipiplet is maintained throughout,
%      irrespective whether lower discarded global eigenstates exist
%      in other symmetry sectors also targeted.
%
% Wb,Apr08,14 ; Wb,Jul06,16

% added possibility to keep NPsi low-energy states simultaneously
% for former version, see: Archive/initNRG_160728.m
% Wb,Jul06,16 // Kyoto

% TST clear; wsys='SU2_w2'; NPsi=1; dbwrn; tst_Hamilton1D
% TST clear; wsys='SU2_w2'; NPsi=4; Qtot=[]; dbwrn; tstflag=2; tst_Hamilton1D
  Iout=struct; tuneH=0; ot={};

  [rs,nq]=symrank(HAM.oez(1).op);
  ns=numel(rs);

  getopt('INIT',varargin);
     Nkeep=getopt('Nkeep',16);
     maxS =getopt('--maxS');
     cplx =getopt('--cplx');

     dQtotN=getopt('dQtotN',[]);
     if ~isempty(dQtotN), [m,n]=size(dQtotN);
        if n~=nq+1 || size(uniquerows(dQtotN(:,1:end-1)),1)<m
           wberr('invalid dQtotN=[%s]',...
           mat2str2(dQtotN,'fmt','%g','rowsep','; ','-f'));
        end
        NPsi=sum(dQtotN(:,end)); Qtot=[];
     else
        NPsi=getopt('NPsi',0);
        if NPsi>0
             Qtot=[];
        else Qtot=0;
        end
     end

     Qtot=getopt('Qtot',Qtot);
     if ~isempty(Qtot)
        if iscell(Qtot), ot=Qtot(2:end); Qtot=Qtot{1}; tuneH=2;
        elseif getopt('-f'), tuneH=2;
        else
           Qopl=getopt('Qop',[]);
           if ~isempty(Qopl), Qop=QSpace; tuneH=1;
           end
        end
     end

     rtol=getopt('rtol',1E-8);

     if getopt('-v'), vflag=1;
     elseif getopt('-V'), vflag=2; else vflag=0; end
  getopt('check_error');

  gotsym=~isempty(get_sym_info(HAM));

  consistency_check(HAM);
  addto_dmrg_info(HAM,HAM);

  L=numel(HAM.mpo);
  E=get_ops_E(HAM,1,'-t'); nq=size(E.Q{1},2);

  if ~isempty(dQtotN)
     if ~isnumeric(dQtotN) || size(dQtotN,2)~=nq+1 || any(dQtotN(:,end)<=0)
        dQtotN, wberr('got invalid dQtotN (%g/%g)',size(dQtotN,2),nq);
     end
     s=sprintf('dQtotN=[%s] @ n={%s}', ...
       mat2str2(dQtotN(:,1:end-1),'fmt','%g','rowsep','; ','-f'), ...
       vec2str(dQtotN(:,end),'sep',',','-f') ...
     );
  elseif ~isempty(Qtot)
     if isequal(Qtot,0), Qtot=zeros(1,size(E.Q{1},2));
     elseif ~isnumeric(Qtot) || ~isequal(size(Qtot),[1, nq]), Qtot
        wberr('got invalid Qtot (%g/%g)',numel(Qtot),nq);
     end
     s=[ 'Qtot=[' vec2str(Qtot,'-f') ']' ];
  else s='Qtot=[]';
  end

  if tuneH, s=[s sprintf(' w/tuneH=%g',tuneH)]; end

  wblog('---',repmat('-',1,50));
  if NPsi>1
       wblog('<i>','NRG initalization of %g DMRG (ground) states',NPsi);
  else wblog('<i>','NRG initalization of DMRG ground state');
  end
  if length(s)>20
  wblog(' * ','%s',s); s=''; else s=[s ', ']; end
  wblog(' * ','%sNkeep=%g, rtol=%g',s,Nkeep,rtol);
  wblog('---',repmat('-',1,50));

  found=0;

  for k=L:-1:1
     E=get_ops_E(HAM,k,'-t');
     Xk=load_dmrg_data(HAM,k);

     if k<L
        Xk.AK=QSpace(permuteQS(getIdentityQS(Xl.AK,1,E),[3 1 2]));
     else
        Xk.AK=getAtensorLoc(E,[],'R');
        q=use_Hconj(HAM);
        if q, q(2)=Xk.info.hconj;
           if xor(q(1),q(2)), wberr('inconsistent hconj setting (%/%g)',q);
           elseif any(q<0), wberr('invalid hconj setting (%/%g)',q); end
        end
     end

     Xk.AK=setitags(Xk.AK,'-A',k);

     if vflag>1
        d=getDimQS(Xk.AK);
        wblog(' * ','AK: [%s] @ [%s]',vec2str(d(1,:)), vec2str(d(end,:))); 
     end

     if ~gotsym && tuneH==1, t3=sprintf('s%02g',k); Qop_=Qop;
        Qop=QSpace(contractQS(Xk.AK,'!1*',{Xk.AK,Qop})) ...
          + contractQS(Xk.AK,'!1*',{Xk.AK,Qopl,['-op:' t3]});
     end

     if k<L
        Xk=updateHK(HAM,k,'<<',Xk);

        HKt=getBlockHK(Xk,'--fix'); % HK `t'uned

        if k>1, D=Nkeep; else D=-1; Iout.H0=HKt; end

        if tuneH && k>1
           if tuneH==1
              Qkt=Qtot*((L-k)/(L-1));
              ot={Qop};
           else Qkt=Qtot; end

           [HKt,It(k)]=tuneHamQtot(HKt,Qkt,ot{:});
        end

     end
     if k<L && HKt

      % NB! HK is rebuilt below via HKt -> Xk.AK using updateHK()
      % with EKt used for logging purposes and eventually for k==1 but
      % mostly ignored otherwise [together with HKt -> tuneHamQtot()!]
      % --> no need to undo tuning of Hamiltonian, since eventually
      % only HAM is used to build HK! --> can ignore UNDO_TUNE_QOP
      % Wb,Feb02,23
        [ee,Ie]=eigQS(HKt,'Nkeep',D); EKt=Ie.EK;
        Xk.AK=contract(Ie.AK,1,Xk.AK,1);

        nd=numel(EKt.data);
           q0=zeros(1,nd);
           for i=1:nd, q0=min(EKt.data{i}); end
        q0=EKt.Q{1}(find(q0<=min(q0)+1E-12),:);

        qk=EKt.Q{1}; qk=mat2cell(qk,size(qk,1),rs);
        qd=zeros(1,ns);

        for i=1:ns
           qk{i}=uniquerows(qk{i});
           qd(i)=size(qk{i},1);
        end

        if vflag || any(qd<=1), d=getDimQS(Ie.AK);

           if vflag>1, f={'%3g','%4g'}; else f={'%2g','%3g'}; end
           s=sprintf([f{1} '/' f{1}],d(1,[2 1]));
           if size(d,1)>1
              s=[s ' ' sprintf(f{2},d(end,1))];
           end

           s=sprintf('%5g  %+-10.6g  %12s %9s   %s', k, ee(1)/(L-k+1), ...
            ['(' mat2str2(q0,'fmt','%g','rowsep','; ','-f') ')'],...
            ['{' vec2str(qd,'-f','sep',',') '}'], s);

           if gotsym && any(qd<=1)
              s=[s '  ']; if D>2, s=[s 'WRN ']; end
              s=[s 'single Q sector'];
              if tuneH && k<L, s=[s sprintf(' (tuneH=%g)',It(k).wtune)]; end
           elseif tuneH && k<L, q=It(k).wtune;
              if ~isempty(q) && q>0
                 if isfield(It(k),'Ufac') && ~isempty(It(k).Ufac)
                      s=[s sprintf('  tuneH @ U=%.4g',It(k).Ufac)];
                 else s=[s sprintf('  tuneH=%g',q)];
                 end
              end
           end

           if vflag>1
              wblog('%s',s); % if k<=1, fprintf(1,'\n'); end
           else
              if k==L-1, fprintf(1,...
                '\n  nrg_k  energy_e0     qset_q0 w/num-qsectors DK/DX (DK*)\n\n');
              end
              fprintf(1,'  %s\n',s);
           end
        end

        if vflag>1
           wblog(' * ','AK: [%s] @ [%s]',vec2str(d(1,:),'-f'), ...
           vec2str(d(end,:),'-f')); 
        end
     end

     Xk.AK=setitags(Xk.AK,'-A',k);

     if ~gotsym && tuneH==1
        Qop=QSpace(contractQS(Xk.AK,'!1*',{Xk.AK,Qop_})) ...
          + contractQS(Xk.AK,'!1*',{Xk.AK,Qopl,['-op:' t3]});
     end

     if k==1,
        NPsi1=max(1,NPsi); HK=getBlockHK(Xk);
        ltag=iff(NPsi>0,'PSI','WRN');

        eN=ee(:,1); q=1E-6*mean(diff(eN(1:min(end,10))));
        i=find(abs(eN-min(eN))<q); g=numel(i);

        if ~isempty(dQtotN)
           if isempty(Qtot)
              nd=numel(EKt.data);
              e0=zeros(1,nd); for i=1:nd, e0(i)=min(EKt.data{i}); end
              i=find(e0<=min(e0)+1E-12);
              if numel(i)>1, wblog('WRN',...
                'got %g degenerate ground state symmetry sectors',numel(i));
              end
              Qtot=EKt.Q{1}(i(1),:);
           elseif tuneH
              Qtot=round(Qtot);
           end

           QTOT=repmat(Qtot,size(dQtotN,1),1)+dQtotN(:,1:end-1);

           [i1,i2,Im]=matchIndex(EKt.Q{1},QTOT);
           n=[numel(i1), size(dQtotN,1)];
           if n(1)
              if diff(n), wblog('WRN',...
                 'only %g/%g symmetry sectors found for dQtotN @ k=%d',n,k);
              end

              EKt=getsub(QSpace(EKt),i1); ee=sort([EKt.data{:}]);
              for i=1:numel(i1)
                  EKt.data{i}=EKt.data{i}(1:min(end,dQtotN(i2(i),end)));
              end
              q=EKt.data;
              NPsi=numel([q{:}]); NPsi1=NPsi; found=2;

           else Xk.AK, wblog('ERR',...
             'invalid dQtotN=[%s] for NRG low-energy sector (L=%g)', ...
              mat2str2(dQtotN,'fmt','%g','rowsep','; ','-f'),L);
           end
        elseif ~isempty(Qtot)
           if tuneH, Qtot=round(Qtot); end
           i=find(~isnan(Qtot));
           [i1,i2,Im]=matchIndex(EKt.Q{1}(:,i),Qtot(i)); n=numel(i1);
           if n==1
              EKt=getsub(QSpace(EKt),i1); ee=EKt.data{1}; found=1;
              EKt.data{1}=ee(1:NPsi1);
           elseif numel(i)<numel(Qtot)
              i1=matchIndex(HK.Q{1}(:,i),Qtot(i));
              HK=getsub(HK,i1);
           else s=vec2str(Qtot,'-f');
              if n>1, wberr('Qtot=[%s] found %g times in H0 !?',s,n);
              else Xk.AK, wblog(...
                'ERR','Qtot=[%s] not found in H0 (L=%g)',s,L);
              end
           end
        end

        if ~found
           [ee,Ie]=eigQS(HK,'Nkeep',NPsi1); ee=ee(:,1);
           EKt=Ie.EK;

           Xk.AK=contract(Ie.AK,1,Xk.AK,1);
           Xk.EK=EKt;

           if NPsi>0
              display(QSpace(EKt),' ');
           end

           q=EKt.Q{1}; l=size(q,1);
           if l==1, wblog('==>',...
             'got Qtot=[%s] for NRG ground state sector (L=%g)', ...
              vec2str(q,'-f'),L);
           elseif l>1, wblog(ltag,...
             'keeping %g low-energy symmetry sectors %s',l,mat2str(q));
           end
        end

        d=getDimQS(EKt); d=d(:,2);
        if d(1)~=NPsi1
           wberr('got %g/%g low-energy states !?',d(1),NPsi1); 
        end

        if g>NPsi1, wblog('WRN',...
           'got %g-fold degenerate ground state space (@ %.3g)',...
           g, norm(diff(eN(1:g)))/max(1,norm(eN(1:g))));
        elseif d(1)>1
           wblog(ltag,'keeping %g global multiplets %s (%g states)',...
           d(1),mat2str(EKt.Q{1}),d(end));
        elseif d(end)>1
           wblog(' * ','got unique ground state multiplet [%s](%g)',...
           vec2str(EKt.Q{1},'-f'),d(end));
        else
           wblog(' * ','got unique ground state [%s] (d=%g)',...
           vec2str(EKt.Q{1},'-f'),d(end));
        end

        [i1,i2,Im]=matchIndex(Xk.AK.Q{1},EKt.Q{1});
        if isempty(i1), EKt, Xk.AK
           wberr('got symmetry sector mismatch !?');
        end

        if g>1 && NPsi1<=1 && maxS
           wblog('NB!','using symmetrized ground state space (--maxS)');
           u=repmat(1/sqrt(g),1,g);
           for i=1:numel(i1), j=i1(i);
              Xk.AK.data{j}=contract(u,Xk.AK.data{j}(1:g,:,:),2,1);
           end
        else
           for i=1:numel(i1)
              j=i1(i); l=numel(EKt.data{i2(i)});
              Xk.AK.data{j}=Xk.AK.data{j}(1:l,:,:);
           end
        end

      % WRN! still, by truncating to lowest NPsi1 multiplets,
      % this might forever truncate local state space for the first site!
      % => keep full i1 here [not just i1(1)!]  // Wb,Mar26,15
      % => use full local state space in HAM.oez for DMRG iteration
      %    see update_psi_2site => get_local_id // Wb,May19,17
        Xk.AK=getsub(QSpace(Xk.AK),i1);

        q=getIdentity(Xk.AK,3,'');
        if ~sameas(q,E,'-l')
            if NPsi1>1
                 s=sprintf('%g multiplets',NPsi1);
            else s='gs-multiplet'; end
            wblog('WRN','using %s reduced local state space!',s);
            display(E,'-c'); fprintf(1,'-->'); display(q,'-c');
        end

        e=QSpace(contractQS(Xk.AK,'23*',Xk.AK,'23'));
        if ~isIdentityQS(e), e
           wberr('got unnormalized A(1) tensor !?');
        end
     end

     if cplx, Ak=Xk.AK.data;
        for i=1:numel(Ak)
            Ak{i}=Ak{i}+1E-14i*randn(size(Ak{i}));
        end
        Xk.AK.data=Ak;
     end

     Xk=updateHK(HAM,k,'<<',Xk);

     save_dmrg_data(HAM,Xk,k);
     if k>1, Xl=Xk; 
        if k<L && HKt, Xl.EK=EKt; end
     end
  end

  H0=getBlockHK(Xk);

  d=getDimQS(Xk.AK); d=d(:,1);

  if d(end)>1 || NPsi>0
     wblog(ltag,'keeping global Psi index with current site');
     Q=getIdentity(getvac(Xk.AK),Xk.AK,1);
     Q=setitags(Q,{'K00','','Psi'});
     Xk.AK=contract(Q,'2',Xk.AK,1,[1 3 4 2]);
  else
     wblog('Psi','got single scalar state -> skip qlabel using 1j');
     q=getIdentity(Xk.AK,1,'-0');
     Xk.AK=contract(q,1,Xk.AK,1);
  end

  k2=min(5,floor(L/2));

  q=Xl; Xl=Xk; Xk=q;
  for l=1:k2, if l>1
     Xk=load_dmrg_data(HAM,l+1); end
     [Xl.AK,Xk.AK,r1,I1]=ortho2site(Xl.AK,Xk.AK,'>>',rtol);
     Xl=updateHK(HAM,l,'>>',Xl);
     save_dmrg_data(HAM,Xl,l);
     if l<k2, Xl=Xk; end
  end

  q=Xl; Xl=Xk; Xk=q;
  for k=k2:-1:1, if k<k2
     Xk=load_dmrg_data(HAM,k); end
     [Xk.AK,Xl.AK,r2,I2]=ortho2site(Xk.AK,Xl.AK,'<<',rtol);
     Xl=updateHK(HAM,k+1,'<<',Xl);
     save_dmrg_data(HAM,Xl,k+1);
     if k>1, Xl=Xk; end
  end

  Xk.info.itags=Xk.AK.info.itags;

  save_dmrg_data(HAM,Xk,1);

  addto_dmrg_info(HAM,NPsi);

  if nargout>1
     Iout=add2struct(Iout,HKt,EKt,Qtot,dQtotN,found,Nkeep,maxS,rtol,Xk,NPsi);
     if tuneH, Iout=add2struct(Iout,'Qop?','Qopl?',It); end
  end

end

% -------------------------------------------------------------------- %

