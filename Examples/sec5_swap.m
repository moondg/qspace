% -------------------------------------------------------------------- %
% Construct swap operator for two spins S=1                            %
% [Sec. 5.2 in documentation]                                          %
% Wb,May10,24                                                          %
% -------------------------------------------------------------------- %

  echo_c % 2S=2 => S=1
  [S,IS]=getLocalSpace('Spin','SU2',2);

  A = getIdentity(S,S,[1 3 2]);  % may also use S -> IS.E instead
  X = contract(A,'13*',A,'31');  % swap operator

  echo_c, A = setitags(A,'-A',2)    % setitags(A,{'K01','K02','s02'})  % same
  echo_c, X = setitags(X,'-op:K',2) % setitags(X,{'K02','K02'})  % same

  SS0 = getIdentity(A,2);
  SAS = contract(S,'-op:K01','*',{A,S,'-op:^s'});
  echo_c  % matrix elements for (S'.S)^1
  SS1 = contract(A,'!2*',SAS)
  echo_c  % matrix elements for (S'.S)^2
  SS2 = contract(A,'!2*',{S,'-op:K01','*',{SAS,S,'-op:^s'}})

% all QSpace tensors here have their records sorted the same way
% i.e., according to symmetry labels => can simply catenate RMTs
% WRN! cannot use x=[X.data{:}] here
% [see discussion with subsref() in QSpace documentation, p. 113]
  x=X.data; x=[x{:}]';
  M=zeros(3); M(:,1)=abs(x); % same as SS0.data
  s=SS1.data; M(:,2)=[s{:}];
  s=SS2.data; M(:,3)=[s{:}];

  echo_c % expansion coefficients of X in terms of (S'.S)^n with n=0,1,2
  a=M\x

% optional: safe conversion to int
% a=chopd(a);  % chop off double precision noise

  echo_c % counter check: same within numerical noise, indeed
  e=norm( a(1)*SS0 + a(2)*SS1 + a(3)*SS2 - X )

