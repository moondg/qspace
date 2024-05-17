% SIAM impurity parameters (all energies in units of the half bandwidth D:=1)
%
%   U        onsite interaction
%   epsd     impurity level position
%   Gamma    hybridization strength of impurity
%   B        magnetic field acting on impurity (Zeeman splitting; applied as -B Sz)

  setdef('U',0.12,'Gamma',0.01,'B',0);  % set default values
  setdef('epsd',-U/2);                  % half-filling by default

% NRG discretization parameters
%
%   Lambda   coarse graining strength (dimensionless; typically >= 2)
%   N        length of Wilson chain
%   Nkeep    number of multiplets to keep (or states if all-abelian)
%   Etrunc   energy truncation threshold (complementary to Nkeep)
%   z        shift of logarithmic discretization (value in [0,1[)
%
% By setting Nkeep large, by default, preference is given to
% truncation by energy; using Nkeep, nevertheless, as safeguard.

  setdef('Lambda',2,'N',55,'Nkeep',1024,'Etrunc',7,'z',0);

% Collect model and NRG parameters (global info structure for reference)
  global param
  param=add2struct('-',U,epsd,Gamma,B,N,Lambda,z);

% NRG hopping parameters
  ff=getNRGcoupling(Gamma,Lambda,N,'z',z); % default z=0 if not specified

% Local state space: single spinfull fermionic level
% For the SIAM this is the same for the impurity as well as the Wilson sites.
% The option '-v' enables more verbose output.

  if ~isset('B')
       SYM='Acharge,SU2spin';
  else SYM='Acharge,Aspin';  % results in SS = [S+, S−, Sz] below => SS(end)=Sz
  end

  [FF,Z,SS,IS]=getLocalSpace('FermionS',SYM,'-v');

% Construct local occupation operator
% depending on SYM above, FF may contain different number of operators (nF)
  N0=QSpace(size(FF)); nF=numel(FF);
     for i=1:nF                                    % sum over i=spin for nF>1
     N0(i)=contract(FF(i),'13*',FF(i),'13'); end   % Fi'*Fi
  N0=sum(N0);

% A-tensor for H0 (LRs index order convention)
  A0 = getIdentity(getvac(IS.E),IS.E,[1 3 2]);

% Impurity Hamiltonian [in R-basis of A0, as in LRs]
  H0 = epsd*N0 + (U/2)*N0*(N0-IS.E) + 0*IS.E;
  if isset('B'), H0 = H0 - B*SS(end); end  % SS(3=end) contains Sz

% Rescale hoppings from global to iterative units of order 1 (see usage below)
  rL=sqrt(Lambda);  % rL = `root Lambda'
  ff=ff.* rL.^(1:N-1);

% Save data by default to files in local Matlab data directory (LMA)
% as $LMA/NRG/myNRG_##.mat with ## the Wilson shell index k
% where LMA is expected to exist as environmental variable
  setdef('sflag',1);  % set sflag=0 to disable
  if sflag
     fout={'NRG' '/myNRG'};
     cto lma;  % `change to' directory $LMA
     if ~isdir(['./' fout{1}]), system(['mkdir ' fout{1}]); end
     fout=[fout{:}];
     wblog('I/O','using %s_##.mat',fout); 
  end

% Initialize iterative diagonalization
%    K  kept space      e.g., AK
%    D  discarded space e.g., AD
% scalar operators only have KK -> K, DD -> D (like HK, HD)
% general operators also have off-diagonal blocks KD, DK (like FKK)

  k=0; % start at Wilson shell `0' (impurity only)
  AK=setitags(setitags(A0,'-A',k),1,'Lvac');   % sets itags {'Lvac','K00','s00'}
  HK=setitags(H0,'-op:K',k);                   % sets itags {'K00','K00*'}
  AD=QSpace;                                   % inits to empty QSpace [= QSpace()]

% Fermionic operator F in KK space
  FKK=QSpace(1,nF);

% Keep track of general data along iterations
  EX=nan(N,2);     % 2 columns: heighest kept and lowest discarded energy
  E0=zeros(1,N);   % subtracted ground state energy for each iteration
  N4=[];           % number of kept and discarded multiplets/states

% Initialize structure to collect iterative diagonalization
  Inrg=struct('N',N,'Lambda',Lambda,'EK',[],'HK',QSpace(1,N));

  fprintf(1,'\n'); disp(param);
  fprintf(1,'\n');

% Iterative diagonalization based on couplings ff(k)
  for k=1:N  % Wilson shell index
     if k<N || ~sflag
          o={'Nkeep',Nkeep,'Etrunc',Etrunc};
     else o={'Nkeep',0};  % by default, discard all states at last iteration
     end

   % Collect info  on total expanded dimension (K+D)
   % this copies 2 entries for non-abelian: [multiplet, state space] dimension
   % where q(1,*) are multiplet dimensions
     q=getDimQS(HK); N4(k,:,2)=q(:,2); % total dimension K+D

   % Exact diagonalization of expanded space (may also use QSpace/eig wrapper here)
     [ee,I]=eigQS(HK,o{:});
     ee=ee(:,1);         % first column (relevant for non-abelian symmetries only)
     EK=QSpace(I.EK);    % MEX files cannot return QSpace objects, only structures
     ED=QSpace(I.ED);    % hence the conversion to QSpace objects here

   % Subtract ground state energy [except for k==1: H0 -> keep E0(1)=0]
     if k>1, E0(k)=min(ee); ee=ee-E0(k); end

     if EK % i.e., non-empty
      % EK is returned in compact diagonal format => take dimension on 2nd index
      % Dk(1,*) is multiplet dimension
        Dk=getDimQS(EK); N4(k,:,1)=Dk(:,2);
        EK=EK-E0(k); q=EK.data;  % also subtract E0 in QSpace EK
        EX(k,1)=max([q{:}]);     % largest kept energy
     end
     if ED % repeat for ED
        ED=ED-E0(k); q=ED.data;
        EX(k,2)=min([q{:}]);     % smallest discarded energy
     end
     
     if sflag % obtain AD before AK is changed right next
     AD=contract(AK,I.AD,[1 3 2]); end
     AK=contract(AK,I.AK,[1 3 2]); % LRs index order convention
     if sflag
        q=struct('AK',AK,'AD',AD,'HK',EK,'HD',ED,'E0',E0(k));
        save(sprintf('%s_%02g.mat',fout,k-1),'-struct','q');
     end

   % Collect `finite size' energy spectra
   % (note that nrg_plot.m prefers Inrg.HK if present, plots EE otherwise)
     m=min(Nkeep,length(ee));
     EE(1:m,k)=ee(1:m); % plain energies without symmetry resolution
     Inrg.HK(k)=EK;     % same energies but with symmetry resolution

     if k==N, fprintf(1,'\n\n'); break; end

   % Wilson shall k actually starts here
   % compute matrix elements for next iteration (propagate FKK)
     for i=1:nF
      % operator index order KK[op] here by having FF to the right
        FKK(i)=contract(AK,'!2*',{AK,FF(i),'-op:^s'});
     end

   % Add new shell (described by `local space' in IS.E)
     AK=getIdentity(I.AK,2,IS.E,[1 3 2]); % LRs convention
     AK=setitags(AK,'-A',k); % generates itags { K<k-1>, K<k>*, s<k>}

   % Rescale and propagate Hamiltonian
     HK=diag(EK*rL); % expand compact diagonal representation
     HK=contract(AK,'!2*',{HK,AK});

   % add hopping to newly added site
     for i=1:nF
        q=ff(k)*Z*FF(i);  % include fermionic parity Z here!
        Q=contract(AK,'!2*',{FKK(i),'*',{AK,q,'-op:^s'}});
        HK=HK+Q+Q';  % +Q' adds hermitian conjugate (" + H.c.")
     end

   % Make sure, zero-diagonal blocks are also included (important for k=1 only)
     HK=HK+0*getIdentity(AK,2);

   % Generic log output: current Wilson shell @ largest kept energy
   % (to be compared to chosen truncation energy Etrunc)
     fprintf(1,' %4d/%d (%g) EK=%.3g ... \r',k,N,max(Dk(1,:)),EX(k,1));

  end % of iteration (Wilson shell) k

% Finalize Inrg info structure (used by nrg_plot.m)
  Inrg.E0=E0;
  Inrg.EK=EX;
  Inrg.EScale=Lambda.^(-(0:length(E0)-1)/2);
  Inrg.Itr=add2struct('-',Etrunc,Nkeep);  % Itr = info on truncation

% Cumulative `physical' value for E0
% i.e., global ground state energy in units of bandwidth
  Inrg.phE0=sum(E0.*Inrg.EScale);

  if ndims(N4)==2
       Inrg.NK=N4;   % abelian symmetries only
  else Inrg.NK=reshape(N4,[],4); % non-abelian
  end

  if sflag
     f=sprintf('%s_info.mat',fout);
     save(f,'-struct','Inrg');
  end

% Finalize numeric array EE for finite-size spectra
% collected data sets have variable length => replace trailing zeros
% for shorter records by nan except so that these are ignored by plot()
  EE(~EE)=nan;            % same as EE(find(E==0))=nan;
  EE(1,isnan(EE(1,:)))=0; % restore zero for first value (ground state)

% Summarize result by plotting NRG energy flow diagram
  if ~exist('plotflag','var') || plotflag  % disable by setting plotflag=0
     param.D=Nkeep;
     nrg_qdisp=10;  % turn on legend with symmetry labels (at most 10 entries)
     nrg_plot       % main plot script for NRG energy flow diagram
  end

  return % comment out to continue

  if sflag
   % Compute the dynamical correlation function <A|B'> where A=op1 and B=op2
   % are operators acting on the impurity (more precisely, acting within A0)
     op1=[];    % operator A (empty takes default: op1=op2)
     op2=FF(1); % operator B: here fermionic annihilation operator
     zflags=1;  % whether to use fermionic parity Z with op[12] (here: yes)

     nostore=1; % do not store matrix elements of op[12] in NRG basis to files fout*
     locRho=1;  % use local thermal density matrix, i.e., do not store Rho(T) 
                % to files fout* either, but keep it in memory

   % Look for relevant options set in Matlab for fdmNRG_QS
   % where a trailing '?' indicates a non-mandatory option
     ofdm=setopts('-','T?','zflags?','-nostore?','-locRho?','nlog?','emin?','emax?');
   % main fdm-NRG routine
     [om,a0,Ifdm] = fdmNRG_QS(fout,op1,op2,Z, ofdm{:});
     fdm_plot;    % plot script to summarize its output
  end


