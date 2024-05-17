
% select what to do / show by setting task (string)
  setdef('task','SU2');
%    'SU2'    spin-half with U1 symmetry
%    'U1'     spin-half with U1 symmetry
%    'Z2'     spin-half with Z2 symmetry
%    'nosym'  spin-half with no symmetry
% Wb,May10,24

  nfin=0; % number of tasks finished

switch task
case {'SU2'}

% using spin operator that describe the spin-half
% using no symmetry at all ('nosym')
  echo_c([task ' spin-half']);
  [S,IS] = getLocalSpace('Spin',1/2);
   % S   % spin operator
   % IS  % info structure, e.g., contains identity operator in IS.E

  S

  echo_c, S2=contract(S,'13*',S,'13') % Casimir operator

  nfin=nfin+1;
end % nosym

switch task
case {'U1','Z2'}

% using spin operator that describe the spin-half
  if isequal(task,'U1')
       sym='-A'; % using U(1) symmmetry [abelian 'A' is synonymous with U(1)]
  else sym='--Z2'; % using Z_2 parity symmetry
  end

  echo_c([task ' spin-half']);
  [S,IS] = getLocalSpace('Spin',1/2,sym);
   % S   % spin operator
   % IS  % info structure, e.g., contains identity operator in IS.E

  S

  echo_c, Sop=sum(S) % combine to full spin operator (no longer an irop)
  echo_c, S2=contract(Sop,'13*',Sop,'13') % Casimir operator

  nfin=nfin+1;
end % U1, Z2

switch task
case {'nosym'}

% using spin operator that describe the spin-half
% using no symmetry at all ('nosym')
  echo_c([task ' spin-half']);
  [S,IS] = getLocalSpace('Spin',1/2,'--nosym');
   % S   % spin operator
   % IS  % info structure, e.g., contains identity operator in IS.E

  S

% this is no longer meaning full here
% Sop=sum(S) 
  echo_c, S(1).data{1} % Pauli  sigma_z / 2
  echo_c, S(2).data{1} % Pauli  sigma_- / sqrt(2)
  echo_c, S(3).data{1} % Pauli -sigma_+ / sqrt(2)

  nfin=nfin+1;
end % nosym

if ~nfin
   if ischar(task)
        printf(1,'\n   ERR invalid task ''%s''\n\n',task);
   else printf(1,'\n   ERR invalid task\n\n'); disp(task);
   end
end

