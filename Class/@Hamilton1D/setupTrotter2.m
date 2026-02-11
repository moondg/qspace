function [HAM,Iout]=setupTrotter2(HAM,dt,varargin)
% function [HAM,Iout]=setupTrotter2(HAM,dt [,opts])
%
%   Setup 2nd order Trotter MPO for specified dt,
%   only accepting Hamiltonian with nearest-neighbor (NN) terms.
%
%   This also initialize the system for Trotter2 time evolution
%   (sweepPsi_tmpo), taking given (ground) state in HAM->AK,
%   and applying op0, i.e, 'ops' locally at specified site k0,
%       |A(t=0)> = op0*|A0> / norm, i.e., without dagger,
%   and subsequent normalization (!)
%   e.g., for a simple normalization of latter truncation error etc.
%
%   This is used with HAM/sweepPsi_tmpo.m to compute the real-time
%   evoultion |At> = exp(-i(H-E0)t) * op0*|A0>, together with the
%   structure factor, i.e., correlation function
%   <A0|op_i' * |At> for all sites i.
%
% Options
%
%    'ops',{X1,k1,X2,k2,...}  if Xi is a set, only first operator is taken
%    'ops',{{X1,iop},k1}      if X1 is a set
%              apply specified operator(s) to given state prior to
%              performing time step; while `one' operator is expected,
%              it may effect a cluster of sites e.g. current, chiral, etc.
%
%    '-f'      enforce MPO setup from scratch
%              i.e., restart from ground state (requires `ops')
%    '--store' store Psi(t) in file structure [by default, keep in memory with HAM]
%
% Wb,May31,19

% -------------------------------------------------------------------- %
% Data structure
%
%     The Trotter data has been assigned a separate variable space
%     under HAM.user.trotter in order to facilitate an indepdendent
%     access to the time dependent wave function, without having
%     to e.g. load/store the ground state wave function every time 
%     for every substep when sweeping along the chain.
%
%     HAM.user.trotter.info    data required/generated with setup of ...
%     HAM.user.trotter.mpo  -> see setup_mpo_trotter2(HAM,dt) below
%
%     HAM.user.trotter.mat     pointer to file structure if ntermediate
%                              data is stored on hard drive; else this
%                              is stored with HAM under
%     HAM.user.trotter.data    which is a structure array with fields:
%
%       .A0  contains original (ground-state) wave function Psi := |0>
%            required for expectation values (such as structure factor)
%       .At  is Psi(t) := e^(-iHt) S_center |0>
%       .Af  is the fitted (f) wavefunction for the next time step t -> t+dt
%            once its finalized (f), this becomes the new At
%            [e.g. see acceptTrotter2Step()]
%       .Xm  block-overlap matrices <Af|..|At> required for
%            the fitting process to obtain the new Af!
%
%  Itag label convention: using `Knn' for A0 and At, `Fnn' for Af
%  for the following reason // tags: AF_ITAGS // Wb,Jul02,19
%
%  In the optimizing the overlap Xm between <Af|Ut|At>, At is constant.
%  (only if At carries a non-scalar global (symmetry) label q_global,
%  the orthogonalization center of At must be carried along, in which
%  case At is not exactly constant). If the operator acting on A0
%  has singleton operator dimension (such as `Sz', or `S_+', etc.)
%  then At can be left as it is. In this case, the overlap matrices Xm 
%  may end up with the SAME bond label `Knn' and the same conj-flag
%  (orientation) for the bon-indices! 
%  
%  To avoid this sitution, the itag labelling  convention for Af
%  is changed using K -> F, i.e. rather than {'Knn','Kn+1','sn+1'}
%  as for At, it uses  {'Fnn','Fn+1','sn+1'} (symbolic notation)
%
% -------------------------------------------------------------------- %

  getopt('init',varargin);
   % dt must be present / is not an option
     ops  =getopt('ops',{});
     force=getopt('-f');
     store=getopt('--store');

     if getopt('-v'), vflag=2;
     elseif getopt('-q'), vflag=0; else vflag=1; end

  getopt('check_error');

  L=length(HAM); kops=0; nops=numel(ops)/2;
  fop=0; ops_in=ops; Sx0=[];

  if ~nops, if force
     wbdie('got empty ops starting from ground state !?'); end
  else
     if ~iscell(ops) || mod(nops,1), wbdie('invalid ops'); end
     if nops~=1
        wbdie('not yet implemented/tested for multiple ops');
     end

     for i=1:nops, j=2*i-1;
        if iscell(ops{j}), q=ops{j};
            ops_in{j}=q{1};
            ops{j}=q{1}(q{2});
        end
        q=isFermOp(HAM,ops{j}); if any(q(:)), fop=1;
           if nops~=1 || numel(q)~=1 || numel(ops{2*i})~=1
           wbdie('invalid fermionic operator set'); end
        end
     end
     ops_=ops;

     if fop
        Zop=HAM.oez(2,1).op; k=ops{2}; fop=k; kops=L;
        ops=[QSpace(1,k-1),ops{1},repmat(Zop,1,L-k)];
     else
        kops=max([ops{2:2:end}]); ops=QSpace(1,kops);
        for i=1:2:numel(ops_)
           ops(ops_{i+1})=ops_{i};
        end
     end
  end

  if isempty(HAM.user), HAM.user=struct; end
  if ~isfield(HAM.user,'trotter'), q=QSpace(1,L);
     if store
        if isempty(HAM.mat), wbdie(...
          'invalid usage (Trotter requires DMRG file store)'); end
        x=[]; mat=[HAM.mat '-dt'];
     else x=q; mat=[]; end

     HAM.user.trotter=struct('info',[],'mpo',q,'data',[],'mat',mat);
  end

  q={ getfield2(HAM.user.trotter,'info','dt', {[]})
      getfield2(HAM.user.trotter,'info','ops',{[]}) };

  if ~force && isequal(dt,q{1}) && (isempty(ops_) || isequal(ops_,q{2}))
     if nargout<2, return; end
  end

  wblog('TR2','setup 2nd order Trotter-Suzuki (dt=%g)',dt);
  HAM=setup_mpo_trotter2(HAM,dt);
  kloc=getCurrentSite(HAM);
  if kloc>1, wbdie('got current site at k=%g/%g !?',kloc,L); end
  kmax=max(kloc,kops);

  if fop, wblog('NB!','applying fermionic operator'); end

  if nargout>1
  XL=QSpace(1,L); XR=QSpace(1,L); end

  for k=L:-1:1
     if k==L, Ik=load_dmrg_data(HAM,k  ); At=Ik.AK; else Ik=In; At=An; end
     if k>1,  In=load_dmrg_data(HAM,k-1); An=In.AK; end

     if kmax<k && itags2odir(At)~=1, kmax=k; end
     q=0; s=sprintf(' \r');

     if k<=kops && ~isempty(ops(k))
        At=contract(At,ops(k),'-op:^s');
        if vflag>1 || fop<=0 || k<=fop, q=1; s={'',''};
           if k<=fop, s{1}=' fermionic'; end
           s{2}=isOpConj(ops(k),'-dag');
           s=sprintf('(applying local%s operator %s)',s{:});
        end
     end

     if vflag>1
        fprintf(1,'\r   k=%2g/%g load data / initialize %s\n',k,L,s);
     elseif vflag && q
        wblog('--> k=%2g/%g %s',k,L,regexprep(s,'[()]',''));
     end

   % orthonormalize (if necessary)
   % NB! whether original Psi itself or any of the operators
   % carries an irop index, need to account for possible 
   % global index/indices; NB! if (multiple) ops are present,
   % all Psi/operator indices will be kept trailing;
   % ops are contracted making use of itags, so for multiple
   % operators, the final state may have a single irop index,
   % which is merged with a global index in Psi (if present) at k=1

     if k>1, if k<=kmax
        [An,At]=ortho2site(An,At,'<<','-l'); end
     else
        q=norm(At); e=abs(q-1);
        if vflag>1
           fprintf(1,'\r%60s',''); if kops || fop, fprintf(1,'\n'); end
        end
        if e>1E-14
           wblog(' *  normalizing state (%.4g -> 1)',q);
           At=(1/q)*At;
        else
           wblog('WRN','state remained normalized @ %.3g !?',1-q);
        end

        r1=numel(At.Q);
        if r1>3, i=finditag(At,'^[^sK]'); % 'op'
           if i==r1 && r1==4
              At=setitags(At,'gop',4);
           else wbdie('unexpected itags for rank-% A-tensor at k=1',r1);
           end
        end
     end

     Af=At;

     if k>1, if k<L, l=1:2; else l=1; end
     else    l=2; end
     Af=itagrep(Af,l,'K','F');

     ut=HAM.user.trotter.mpo(k);
     if k<L
          Q={ut,Xm}; if k==1, Xm=[]; end
     else Q=ut; Xm=[];
     end

     A0=Ik.AK;

     HAM=save_trotter_data(HAM, add2struct('-',A0,Af,At,Xm), k);

     Xm=contract(Af,'!1*',{At,Q}); % LL',m(po)

     if nargout>1
        Q=Af; if fop, Q={Q,Zop,'op:^s'}; end
        if k<L, Q={Q,XR(k+1)}; end
        XR(k)=contract(A0,'!1*',Q);
     end

     Il=Ik;
  end

  HAM.user.trotter.info.ops=ops_in;
  HAM.user.trotter.info.fermionic=fop;

  if vflag, z=getscalar(Xm);
     x=real(z); x(2)=round(x); x=[x(2),x(1)-x(2),imag(z)];
     if nops, s='ops'; else s=''; end
     s=sprintf('<Psi| U(dt) %s|Psi> / nrm2 = %.3g %+3g %+.3gi',s,x);
     if abs(z-1)<1E-14, s=[s ' (unchanged!)']; end
     wblog(' * ',s); 
  end

  if nargout<2, return; end

  if nops
     wblog('t=0','computing correlation function Sx0'); 
     Sx0=zeros(1,L);
     for k=1:L
        Il=load_trotter_data(HAM,k);

        Af=Il.Af; if k>1, Af=contract(XL,Af); end
        Qf=Af; if k<L, Qf={Qf,XR(k+1)}; end
        x=contract(Il.A0,'*',{Qf,ops_{1},'*','-op:^s:gop'});
        Sx0(k)=getscalar(x);

        XL=contract(Il.A0,'!2*',Af);
     end
  end

  Iout=add2struct('-',Sx0,ops,dt);
  Iout.normAt=q;

end

