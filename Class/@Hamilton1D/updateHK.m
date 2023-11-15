function [Xk,e0]=updateHK(HAM,kc,kdir,varargin)
% function [Xk,e0]=updateHK(HAM,kc,kdir[,Xk,Xn])
%
%    [Usage #1] Update Hamiltonian (HK) and operator (OP)
%    matrix elements (this loads the data from the previous
%    iteration, in order to propagate the operators correctly).
%
%    If not specified, Xk will be loaded from DMRG_#k#.
%    The updated HK and OPs will be stored (and returned).
%    However, if Xk is explicitly specified, the DMRG data
%    in store will not be changed, rather the updated data
%    for site kc (Xk) is returned as output argument.
%
%    If Xn for the next (current) site is also specified,
%    the reduced density matrix R is computed which allows
%    to return the energy expectation value e0 of the newly
%    added terms to HK.
%
% Usage #2: Xk=updateHK(HAM,kc,kdir,Xk,'-x')  # full mpo only
%
%    Compute partial overlap update as required for controlled
%    bond expansion (CBE). Here honj(i) specifies whether H.c.
%    should be included for Xk.OP(i[,1]) in Xk.Op(i,2).
%
% Wb,Apr08,14

% [08/2023] added usage #2

% ----------------------------------------------------------------- %
% This uses the @Hamilton1D object HAM to update the Hamiltonian
% (either stored in global memory space or on the file system).
% The class object HAM itself is not altered, hence does also
% not get returned => keep it that way! // Wb,Apr10,14
% ----------------------------------------------------------------- %
% Wb,Jan20,16: keep track of expectation value of Hamiltonian terms
% ----------------------------------------------------------------- %
% Note though that the ground state E0:=energy/site in general
% DOES NOT exactly yield sum(e0), i.e.  // Wb,Jan20,16
%
%    E0 =~ E0a:=sum(triu(HSS)) =~ E0b:=sum(tril(HSS))
%
% since the wave function |Psi> is still further optimized by
% continuing to sweep! in general, E0a and E0b are likely to be
% slightly lower than E0, since E0a and E0b contain locally
% optimized term. In contrast, by later further optimization
% of |Psi>, these values are likely to increase, hence leading to
% slightly larger E0. The difference |E0-E0a| or |E0-E0a| should
% be within the current accuracy of the ground state energy.
% ----------------------------------------------------------------- %

   if nargin<3 || nargin>5
      if ~helpthis(nargout,varargin{:}), wbdie('invalid usage'); end
      return
   end

   kdir=check_dir(kdir);

   if using_full_MPO(HAM)
        [Xk,e0]=updateHK_mpo_full  (HAM,kc,kdir,varargin{:});
   else [Xk,e0]=updateHK_mpo_pseudo(HAM,kc,kdir,varargin{:});
   end
end

% -------------------------------------------------------------------- %
function [Xk,e0]=updateHK_mpo_full(HAM,k,kdir,Xk,Xn)

   global gHSS

   loadXk=0; gotR=0; partial=0;
   if nargin<4
      Xk=load_dmrg_data(HAM,k); loadXk=1;
   elseif nargin>4
      if isfield(Xn,'AK'), gotR=1;
      elseif ischar(Xn) && isequal(Xn,'-x'), partial=1;
      else wbdie('invalid argument Xn'); end
   end

   L=numel(HAM.mpo); e0=zeros(0,3); ex=[]; kfix=0; 

   if L<=99, fmt='%02g'; else fmt='%03g'; end
   if kdir>0
        sdir=sprintf(['|' fmt '>'],k); kp=k-1; if kp<0, kp=0; end
   else sdir=sprintf(['<' fmt '|'],k); kp=k+1; if kp>L, kp=0; end
   end

   if kp
      ix=getfield2(HAM.info,'sweep','init_DMRG',{0});
      if numel(ix)==2
         i3=itags2int(Xk.AK);
         if k~=i3(3), wbdie('inconsistent current site k=%g/%g',k,i3(3)); end
         if kdir>0, kp_=abs(i3(1));
         else       kp_=abs(i3(2))+1; end
         Xp=load_dmrg_data(HAM,kp_);

         if kp~=kp_, kfix=1;
         elseif numel(ix)==2 && kp>=ix(1) && kp<=ix(2), kfix=2;
         elseif isempty(Xp.AK) || isempty(Xp.HK), kfix=3; end

         if kfix
            if kfix==1
                if kdir>0 && isequal([kp_+1, kp],ix), kfix=-kfix; kp=kp_; end
            elseif kdir<0
               if kp==ix(1), kfix=-kfix; kp=ix(2)+1;
                  Xp=load_dmrg_data(HAM,kp);
                  if isempty(Xp.AK) || isempty(Xp.HK), kfix=-10*kfix; end
               end
            end
         end
         if kfix>0, s={'',''};
            if kdir>0, s{1}=sprintf('|%g>',k);
            else s{1}=sprintf('<%g|',k); end

            if isempty(Xp.AK), s{2}='empty AK';
            elseif isempty(Xp.HK), s{2}='empty HK';
            else s{2}=sprintf('e=%g',kfix); end

            wbdie('%s unexpected previous site kp=%g / [%g %g] (%s)',...
            s{1},kp,ix,s{2});
         end
      end

      Xp=load_dmrg_data(HAM,kp);
   end

   Ak=QSpace(Xk.AK);

   NPsi=get_NPsi_site(Ak);
   if NPsi && loadXk, wbdie(...
     'invalid usage (must not update HK with rank-4 current site)');
   end

   [q,it]=gotITags(Ak);
   if q<3, it, wbdie('itags required with MPState');
   elseif q>3 && ~NPsi, it, wbdie('got rank-%g A-tensor',q); end
   sloc=it{3}; R=QSpace;

   if gotR
      NPsi=get_NPsi_site(Xn.AK);

      if kdir>0
          if ~NPsi, ic={'!1*'}; else ic={'23*','23'}; end
      elseif ~NPsi, ic={'!2*'}; else ic={'13*','13'}; end

      R=QSpace(contractQS(Xn.AK,ic{1},Xn.AK,ic{2:end}));

   elseif ~partial && all(getqdir(Xk.AK)>0), gotR=-1;
   end

   if isempty(gHSS) && gotR
      gHSS=zeros(L,L);
   end

   nrep=1; oHc={};
   if partial && Xk.info.hconj, nrep=2; end

   for irep=1:nrep
      Q=Ak; if irep>1, oHc={'*'}; end
      if kp, Q={Q,Xp.HK,oHc{:}};
         if irep==1
            if ~Xp.HK, wbdie('got empty Xp.HK'); end
            tp=abs(get_kidx(Xp.HK));
            tk=abs(get_kidx(HAM.mpo(k)));
            if all(tp(3)~=tk(1:2))
               if kfix<0
                  if kdir>0, i=1; else i=2; end
                  Q{2}=setitags(Q{2},HAM.mpo(k).info.itags{i},3);
               else
                 wbdie('got Hamilton MPO inconsistency (%g/[%g %g])',...
                 tp(3), tk(3:4));
               end
            end
         elseif kfix, wbdie('got kfix=%d with partial',kfix); end
         if kdir>0, Q=Q([2:numel(Q), 1]); end
      end
      if kdir>0, ic='!2*'; else ic='!1*'; end

      Q={Q,HAM.mpo(k),oHc{:}}; p=[]; if ~partial, Q={Ak,ic,Q}; end
      Q=contract(Q{:},p);
      if irep==1
           Xk.HK=Q;
      else Xk.HK(irep)=Q; end
   end

   gflag=got_gauge(HAM);
   if gflag
      HG=get_HK_gauge(Ak,Xp,HAM,kdir);
      Q={HG,Ak}; if ~partial, Q={Ak,ic,Q}; end
      Xk.HK = Xk.HK + contract(Q{:});
      if gflag>1, Xk=update_Qop_gauge(Xk,Xp,HAM,kdir); end
   end

   r=numel(Xk.HK(1).Q); l=3; if partial, l=4; end
   if r>l, wbdie('unexpected HK (rank %d/%d)',r,l);
   elseif ~r, wbdie('Hamiltonian MPO contracted to empty',r);
   end

   if gotR
      if gotR>0, q={Xn.info.Eg(end)}; else q={}; end

      HK=getBlockHK(Xk);
      [q,q_]=get_energies(HK,R,q{:});
      if isempty(ex), ex=q_; else ex=ex+q_; end

      e0(end+1,:)=[k,k,q];
      gHSS(k,k,1)=q;

   end

   Xk.info.itags=Xk.AK.info.itags;

   isw=getopt(HAM.info,'sweep?','isw',1);
      if isempty(isw) || isw<1, isw=1; end
   w=get_stype(HAM);

   if gotR
      Xk.info.ek{isw,w}=e0;
      e0=sum(e0(:,end)); if ~gotR && e0==0, e0=nan; end
   end

   Xk.info.ex=ex;

   if loadXk
      save_dmrg_data(HAM,Xk,k);
   elseif ~nargout
      clear Xk
   end

end

% -------------------------------------------------------------------- %
function [Xk,e0]=updateHK_mpo_pseudo(HAM,kc,kdir,Xk,Xn)

   global gHSS

 % NB! no partial contraction here // tags: -x
 % since using a pseudo MPO, this has no MPO bond dimension
 % that can be used; therefore contracting into a 2-site (Ls,Rs)
 % this *is* 2-site => rather used bond expansion within
 % full 2-site setting // Wb,Sep09,23

   gotR=0; loadXk=0;

   if nargin<4, loadXk=1;
   elseif nargin>4
      if isfield(Xn,'AK')
         if nargout>1, gotR=1; end
      else wbdie('invalid argument Xn'); end
   end

   L=numel(HAM.mpo);
   e0=zeros(0,3); ex=[];

   if kdir>0, kp=kc-1; sdir=sprintf('|%02g>',kc);
   else       kp=kc+1; sdir=sprintf('<%02g|',kc); end

   w=HAM.mpo(kc).stype; ioff=size(HAM.oez,1);

     oez = [ HAM.oez(:,w).op        ];
     ops = [ HAM.ops(:,w).op        ];
     dop = [ HAM.ops(:,w).dop       ];
     fop = [ HAM.ops(:,w).fermionic ];
   hconj = [ HAM.ops(:,w).hconj     ];

   M=HAM.mpo(kc).M;
   iop=HAM.mpo(kc).iop;

   s=size(M); s(end+1:3)=1; d=s(3);
   if d<2 || length(s)>3
      wbdie('invalid mpo (got changing d_op=%g @ kc=%g)', d,kc); 
   end

   q=zeros(2,2,d); q(:,:,1)=eye(2);
   x=M(1:2,1:2,:); x(1,2,2:end)=0;
   if ~isequal(q,x), wbdie(...
      'invalid MPO (expecting |0> and |e> to go with first two indices)');
   end

   Q=QSpace;
   if kdir>0 && kc>1 || kdir<0 && kc<L
        Xp=load_dmrg_data(HAM,kp); Q=Xp.HK;
   else Xp=[]; end

   if loadXk
      Xk=load_dmrg_data(HAM,kc );
   end

   Ak=QSpace(Xk.AK);

   NPsi=get_NPsi_site(Ak);

   if NPsi && loadXk, wbdie(...
     'invalid usage (must not update HK with rank-4 current site)');
   end

   [q,it]=gotITags(Ak);
   if q<3, it, wbdie('itags required with MPState');
   elseif q>3 && ~NPsi, it, wbdie('got rank-%g A-tensor !?',q);
   end
   sloc=it{3};

   if gotR
      NPsi=get_NPsi_site(Xn.AK);
      if kdir>0
         if ~NPsi
              R=QSpace(contractQS(Xn.AK,'!1*',Xn.AK));
         else R=QSpace(contractQS(Xn.AK,'23*',Xn.AK,'23'));
         end
      else
         if ~NPsi
              R=QSpace(contractQS(Xn.AK,'!2*',Xn.AK));
         else R=QSpace(contractQS(Xn.AK,'13*',Xn.AK,'13'));
         end
      end
   elseif all(getqdir(Xk.AK)>0)
        gotR=-1; R=QSpace;
   else gotR=0; end

   if isempty(gHSS) && gotR
      gHSS=zeros(L,L);
   end

   if ~isempty(Q), Q={Ak,Q};
        H=contract(Ak,'*',Q,'!1');
   else H=QSpace; end

   gflag=got_gauge(HAM); % isfield(param,'gauge') // Wb,Apr30,23

   if gflag % isfield(param,'gauge') // Wb,Apr30,23
      HG=get_HK_gauge(Ak,Xp,HAM,kdir);
      H=H+contract(Ak,'*',{Ak,HG},'!1');
      if gflag>1, Xk=update_Qop_gauge(Xk,Xp,HAM,kdir); end
   end

 % add local term(s)
 % --------------------------------------------------------------------
 % note that there may be several terms if local operator
 % needs to be written as superposition of other operators
 % such as \hat{n} = (\hat{1} - \hat{z})/2
   I=find(M(1,2,:)); tloc=M(1,2,I);
   if numel(I), I=I-ioff;
      if any(dop(I)~=1)
         wbdie('invalid usage (local ops must be all scalars)'); 
      end

      q=QSpace;
      for i=1:numel(I)
       % allow non-hermitian local operators, but which need to add hconj
       % e.g. ops=[Sz,+S_/sqrt(2)] while skipping S+ = (S_)'
       % ==> Bx*Sx = (Bx/sqrt(2))*ops(2) + H.c.
       % (Sx itself cannot be represented, since Hamilton1D insists
       % on having orthogonal set of local operators HAM.ops // Wb,Aug02,19
         h=tloc(i)*ops(I(i)); if hconj(I(i)), h=h+h'; end
         q=q+h;
      end

      if kdir>0, ic='!2*'; else ic='!1*'; end
      Q=contract(Ak,ic,{Ak,q,['-op:' sloc]});

      if isempty(Q.data)
         wblog('WRN','adding empty local contribution to H');
      end

      H=H+Q;

      if gotR
         if gotR>0, q={Xn.info.Eg(end)}; else q={}; end
         [q,q_]=get_energies(Q,R,q{:});
         if isempty(ex), ex=q_; else ex=ex+q_; end

         e0(end+1,:)=[kc,kc,q];
         gHSS(kc,kc,1)=q;
      end
   end

 % add finalized (completed) coupling terms
 % --------------------------------------------------------------------
 % look for terms that terminate at given site described by M.
 % NB! M matrices are built L->R, moving from |0> to |e>
 % when moving from R->L, hence |e> to |0>!
 % --------------------------------------------------------------------
 % NB! By construction, LR is dealt with *differently* than RL!
 % i.e. consider coupling terms that start at fixed k1 and terminate 
 % at several k2 values (having k1<k2). Then:
 % LR: starting operator at k1 is computed and carried through to 
 % later iterations until it hits any k2, where it is contracted
 % with the corresponding terminating operator; on the other hand,
 % RL: terminating operators are *added* up, and finally 
 % added to H in a *single* contraction at k1 // LR_RL_BEHAVIOR_IN_M
 % --------------------------------------------------------------------

   if kdir>0
        w=1; ic='!2*'; [I,K,t]=find(permute(M(:,2,:),[1 3 2]));
   else w=2; ic='!1*'; [I,K,t]=find(permute(M(1,:,:),[2 3 1]));
   end

   i=find(I>2);
      I=I(i)-2;
      K=K(i)-ioff;
      t=t(i);

   for i=1:numel(I)
      q=t(i)*ops(K(i));
      if fop(K(i)) && kdir>0
         q=oez(2)*q;
      end

      % WORST possible contraction pattern wrt. CData! // Wb,Aug10,16
      % contractQS({ Ak,Xp.OP(I(i)) },'!1*', { Ak,q,['-op:' sloc] });
      % => generates (large) rank-4 objects *twice* *before*
      %    contracting them over 3 indices simultaneously!
      % however, similar contraction cannot be avoided in the case
      % of longer-range interactions in the Hamiltonian (i.e. see
      % update of OP -- if required -- below; WRN_CGC_COST_RANK4)
      Q=contract(Ak,ic,{ Xp.OP(I(i)),'*', { Ak,q,['-op:' sloc] }});

      if hconj(K(i)), Q=Q+Q'; end
      if ~Q, wblog('WRN',...
        'got empty coupling contribution to H (%s)',sdir); end
      H=H+Q;

      if gotR
         if gotR>0, q={Xn.info.Eg(end)}; else q={}; end
         [q,q_]=get_energies(Q,R,q{:});

         if isempty(ex), ex=q_; else ex=ex+q_; end

         [k1,io1]=get_kiop(iop(I(i)+2,w,1));
         [k2,io2]=get_kiop(iop(I(i)+2,w,2));

         if kdir>0
            e0(end+1,:)=[k1,k2,q];
            gHSS(k1,k2,io1)=q;

         else
          % NB! back-sweep is already somewhat convoluted, since
          % contributions are already partially summed here!
          % namely all k2>=k1 to the same k1, with the contribution
          % added here at min(k2) // LR_RL_BEHAVIOR_IN_M
          % => store in lower triangular part of gHSS
            e0(end+1,:)=[k2,k1,q];
            gHSS(k2,k1,io1)=q;

         end
      end

    % ----------------------------------------------------------
    % check whether to include higher powers of coupling term
    % ----------------------------------------------------------
    % NB! powers cannot be simply generated by taking powers
    % of Q above, since the powers in the operators in the
    % representation of the local state space do not necessarily
    % directly carry over to the operator representation
    % in the effective basis, specifically so in the presence
    % of truncation, except for nearest neighbor terms only
    % in the bond picture. // Wb,Apr21,14
    % hcpow is depricated; see QOP_ALLOW_MULT // Wb,Aug26,15
    % ----------------------------------------------------------
   end

   if ~isIdentityCG(H) || ~isHConjQS(H)
      save2('-f','./tmp-HHc.mat'); pwd, e=normQS(H-H'); wbdie(...
     'got non-scalar/non-hermitian Hamiltonian (e=%g @ kc=%g)',e,kc);
   end

   Xk.HK=H;

 % --------------------------------------------------------------------
   if kdir>0
        [J,K,t]=find(permute(M(1,3:end,:),[2 3 1]));
   else [J,K,t]=find(permute(M(3:end,2,:),[1 3 2]));
   end

 % NB! non-unique K cannot be fully avoided!
 % ==> reuse already computed matrix elements
 % Inequivalent coupling terms that start here with the same matrix
 % elements may already be immediatly go different paths by adding
 % (thus merging) certain entries here with propagated terms below!
 % Wb,Aug26,15 // see private/setup_mpo::reduce_RL() // MERGING_TERMS

   [K,i]=sort(K-ioff); J=J(i); t=t(i);

   uJ=all(diff(sort(J(:)))>0);

   if ~isempty(J), n=max(J);
   else n=0; end

   if isfield(Xk,'OP'), n=max(n,numel(Xk.OP));
   else
      n=numel(fieldnames(Xk));
      Xk=orderfields(Xk,[1:n-2, n, n-1]);
   end
   Xk.OP=QSpace(1,n);

 % NB! keep track of originating local operator through itags!?
 % i.e. in addition to mpo(k).iop: ot=sprintf(':op%g',K(i));
 % NB! don't! Heisenberg interaction between sites with
 % different spin or AKLT type Hamiltonian require contraction
 % of *different* ops [see TWIST in setup_Heisenberg()]
 % Wb,Aug27,15

   for i=1:numel(J)
      if i>1 && K(i)==K(i-1) && uJ

         Xk.OP(J(i)) = (t(i)/t(i-1)) * Xk.OP(J(i-1));
      else
         kop=t(i)*ops(K(i));
         if fop(K(i)) && kdir<0
            kop=oez(2)*kop;
         end

         Q=contract(Ak,ic,{Ak,kop,['-op:' sloc]});

         Xk.OP(J(i)) = Xk.OP(J(i)) + Q;
      end
   end

 % --------------------------------------------------------------------
 % propagate operators (if any) and
 % possibly add / merge with newly generated coupling terms above
 % this is required when dealing with longer range-interactions
 % which includes periodic boundary conditions.
 % WRN! this can be costly for larger-rank symmetries
 % since rank-4 CGCs with two large indizes must be generated
 % tags: WRN_CGC_COST_RANK4 // Wb,Aug10,16
 % --------------------------------------------------------------------
   q=M(3:end,3:end,:); I=find(q); t=q(I);
   [I,J,K]=ind2sub(size(q),I);
   if kdir<0
      q=I; I=J; J=q;
   end

   for i=1:numel(I)
      k=K(i); if k>ioff,
           kop=ops(k-ioff);
      else kop=oez(k); end

      if t(i)~=1
         s=sprintf('got mpo weight within operator propagation !? (%.3g)',t(i));
            if k>1, wblog('WRN','%s',s);
            else wbdie('%s',s); end
         kop=t(i)*kop;
      end

      if k>1
           Akop={{Ak,kop,['-op:' sloc]}}; % '*3' enters through op
      else Akop={Ak}; end

      Q={ Akop{:}, Xp.OP(I(i)) };

      try
         Q=contract(Ak,ic,Q);
      catch l
         save2(sprintf('%s_%s_dbg.mat',wbstamp,mfilename));
         save2tmp; rethrow(l)
      end

      if size(Xk.OP)>=J(i) && ~isempty(Xk.OP(J(i)))
           Xk.OP(J(i)) = Xk.OP(J(i)) + Q;
      else Xk.OP(J(i)) = Q;
      end
   end

   Xk.info.itags=Xk.AK.info.itags;

   isw=getopt(HAM.info,'sweep?','isw',1);
   if isempty(isw) || isw<1, isw=1; end

   if gotR
      Xk.info.ek{isw,w}=e0;
      e0=sum(e0(:,end)); if ~gotR && e0==0, e0=nan; end
   end

   Xk.info.ex=ex;

   if loadXk
      save_dmrg_data(HAM,Xk,kc);
   elseif ~nargout
      clear Xk
   end

end

% -------------------------------------------------------------------- %
function [e0,ex]=get_energies(H,R,Eg)

  rH=numel(H.Q);
  rR=numel(R.Q);

  if ~rH || ~rR
     e0=trace(QSpace(H)); if nargout>1, ex=e0; end
  elseif rR==rH
     h=QSpace(contractQS(H,R));
     e0=getscalar(h);     if nargout>1, ex=e0; end
  elseif rR==rH+2
     if isnumeric(Eg)
        q=getIdentityQS(R,4);
        d=getDimQS(q); if d(1)~=1
           wbdie('unexpected R/H setting !?'); end
        R0=getrhoQS(q);
     else
        R0=getrhoQS(Eg);
     end

     h=contract(H,{R,R0});
     e0=getscalar(h);

     if isnumeric(Eg)
        if nargout>1, ex=e0; end
        return
     end

     if nargout>1
        if isdiag(Eg,'-f')==1, Eg=diag(Eg); end
        ee=Eg;
        for k=1:numel(Eg.data)
           q=getsub(Eg,k); n=length(q.data{1}); z=zeros(n);
           [~,~,dc]=getQDimQS(q,1); d=prod(dc);
           for i=1:n
              q.data{1}=z; q.data{1}(i,i)=1/d;

              h=contract(H,{R,q});
              ee.data{k}(i)=getscalar(h);
           end
        end

        if 1
           eg=Eg.data; [eg,is]=sort([eg{:}]);
           ex=ee.data; ex=[ex{:}]; ex=ex(is);
        else
           ex=ee;
        end
     end
  else
     wbdie('unexpected input QSpaces !?');
  end

end

% -------------------------------------------------------------------- %

