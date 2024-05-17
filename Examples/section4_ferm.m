
% select what to do / show by setting task (string)
  setdef('task','overview');
%    'overview'  selected output
%    'struct'    detailed display of data structure of F
%    'noise'     check difference of QSpace tensors w.r.t. numerical noise
%    'Fop'       simple commmands with fermionic annihilation operator F
%    'Sop'       simple commmands with spin operator
%    'perm'      permuting indices of legs and conjugate tensors
%    '1j'        1j tensors
%    '*'         choose all
% Wb,May10,24


% get set of operators that describe the local state space of a
% single spinful site, using U(1) charge and SU(2) spin symmetries
  [F,Z,S,IS]=getLocalSpace('FermionS','Acharge,SU2spin');
   % F   % fermionic annihilation operators
   % Z   % fermionic parity operator
   % S   % spin operator
   % IS  % info structure, e.g., contains identity operator in IS.E

% get identity operator for given local state space
  E=getIdentity(F); % e.g., based on the annihilation operator F here

  nfin=0; % number of tasks finished

switch task
case {'overview'} % this part is run if task='*' or 'struct'
   echo_c('-H','overview fermionic spin-half site');
   echo_c,  F % annihilation operator
   echo_c,  Z % fermionic parity operator
   echo_c,  S % spin operator
   echo_c,  E % identity operator (same as in IS.E)
   echo_c, IS % info structure (mostly for internal purposes, hence ignore)
   nfin=nfin+1;
end

switch task
case {'*','struct'} % this part is run if task='*' or 'struct'
 % explore bare data structure
   echo_c('-h','bare QSpace data structure'); F
   echo_c, struct(F)
   echo_c, F.data'
   echo_c, F.Q{1}, F.Q{2}, F.Q{3}
   echo_c, F.info
   echo_c, F.info.cgr(1,1)
   echo_c, F.info.cgr(1,2)
   nfin=nfin+1;
end

switch task
case {'*','noise'}
   echo_c('-h','adding numerical noise to identity E'), E
   echo_c, E == IS.E  % returns 1 (true), same operators indeed

   echo_c, E.data{1} = E.data{1} + 1E-15; % offset by numerical noise (eps)
   echo_c, E == IS.E  % returns 0 (false) now
   echo_c, norm(E-IS.E) % same output as for (1+1E-15)-1 => 1.1102e-15 (!)
   echo_c
   sameas(E,IS.E)  % allows difference on the order of numerical noise (returns 1=true)

   echo_c, sameas(E,IS.E,1E-16)  % overwrite default 1E-10 => returns 0=false now
   nfin=nfin+1;
end

switch task
case {'*','Fop','nloc'}
   echo_c('-H','simple operations with annihilation operator');
   echo_c, nloc = contract(F,'13*',F,'13') % local occupation operator

   echo_c
   E=acomm(F,F')/2 % anticommutator (using factor 1/2 to account for sum over spin)

   echo_c, isIdentityQS(E) % confirms fermionic anticommutator relation

   echo_c
   [ trace(nloc), norm(F)^2 ] % check norm^2 consistency (expecting same value of 4.)

   echo_c
   X=contract(F,'123*',F,'123') % full contraction with itself also rank-0 tensor
   echo_c, getscalar(X) % that represents the scalar 4. as above

   echo_c
   V=getvac(F,'-1d') % the only possible rank-1 tensor (must have q-labels of vacuum)
   nfin=nfin+1;
end

switch task
case {'*','Sop','spin'}
   echo_c('-H','simple operations with spin operator');
   echo_c, S2=contract(S,'13*',S,'13') % Casimir operator => 3/4 for spin-half
   nfin=nfin+1;
end

switch task
case {'*','perm'}
   echo_c('-H','permutations and index order');

   echo_c
   F2=contract(F,'-op:s01','*', Z*F,'-op:s02', [2 3 1 4])
   % same as      F2=contract(F,'-op:s01','*',Z*F,'-op:s02');
   % followed by  F2=permute(F2,'2314')

   echo_c, F2_ = permute(F2,'3412*') % corresponds to Hermitian conjugate
   echo_c, norm(F2-F2_)  % hence different from F2
   nfin=nfin+1;
end

switch task
case {'*','1j'}
   echo_c('-h','1j tensor and trailing marker characters')

   U=getIdentity(F,'-0')

   echo_c, Uc=conj(U) % conjugate tensor (same index order!)
   echo_c, norm(contract(U,2,Uc,2) - E) % U and Uc are orthogonal matrices indeed
   echo_c, norm(contract(U,2,U,'2*') - E) % equivalent to previous line
   nfin=nfin+1;
end

if ~nfin
   if ischar(task)
        printf(1,'\n   ERR invalid task ''%s''\n\n',task);
   else printf(1,'\n   ERR invalid task\n\n'); disp(task);
   end
end

