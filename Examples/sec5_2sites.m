% -------------------------------------------------------------------- %
% Generate combined state space of two fermionic sites                 %
% [Sec. 5.1 in documentation]                                          %
% Wb,May10,24                                                          %
% -------------------------------------------------------------------- %

% select what to do / show by setting task (string)
%   'intro'   initial setup
%   'SdotS'   compute matrix elements for S.S interaction
%   'FdagF'   Hamiltonian from fermionic hopping term
%   '*'       choose all
  setdef('task','intro');

  nfin=0; % number of tasks finished

  if ~isvar('F' ) || ~isa(F,'QSpace') ||  ...
     ~isvar('IS') || ~isfield(IS,'sym') || ~isequal(IS.sym,'Acharge,SU2spin')

     echo_c('spinfull fermionic site-half');
     [F,Z,S,IS]=getLocalSpace('FermionS','Acharge,SU2spin');

   % fuse state space of 2 sites (LRs index order ^1)
     A=getIdentity(Z,Z,[1 3 2]);

     echo_c
     A=setitags(A,{'K01','K02','s02'})
     % ^1) K01 (kept from `Left', site 1), fused with
     %     s02 (site 2), results in combined
     %     K02 (kept after site 2, `Right')

     echo_c
     E2=getIdentity(A,2)  % identity on combined state space (2 sites)
  end

switch task
case {'*','SdotS'}
% compute S.S interaction
  echo_c
  S12 = contract(A,'!2*',{S,'-op:K01','*',{A,S,'-op:^s'}})

  echo_c, S12.data{1} % singlet space (3 multiplets total)

  nfin=nfin+1;
end

switch task
case {'*','FdagF'}
% compute F'.F interaction (with Fermionic parity Z on second site)
  echo_c
  T12 = contract(A,'!2*',{F,'-op:K01','*',{A,(Z*F),'-op:^s'}});

  echo_c
  H12 = T12 + T12'  % adding H.c. for typical Hamiltonian term

  echo_c % got absent zero-blocks in H12
  skipzeros(E2 - getIdentity(H12),'-f') % skip zero blocks
  % option '-f' enforces skip even for scalar operator

  echo_c % obtain complete eigenspectrum (by adding zero blocks from E2)
  [U12,E12,I12] = eig( H12 + 0*E2 );

  echo_c, E12  % eigenenergies as QSpace (inluces zero-blocks for completeness)

  echo_c, I12.ee  % eigenenergies in numeric format in 1st column

  echo_c  % second column specifies degeneracies
  sum(I12.ee(:,2))  % adds up to 16, as expected

  nfin=nfin+1;
end

if ~nfin
   if ischar(task)
        printf(1,'\n   ERR invalid task ''%s''\n\n',task);
   else printf(1,'\n   ERR invalid task\n\n'); disp(task);
   end
end

