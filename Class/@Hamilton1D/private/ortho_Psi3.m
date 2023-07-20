function [Q,ns,nrm]=ortho_Psi3(Q,stol,NPsi)
% function [Q,ns,nrm]=ortho_Psi3(Q,stol,NPsi)
%
%    Orthonormalize Q in g assuming lr[g] index order.
%    i.e. orthonormalizes states in Q w.r.t to state index (#3)
%    including truncation for states with SVD < stol.
%
%    NB! use SVD-based truncation with state rather than rtol
%    with a density matrix to better separate numerical noise
%    Wb,Feb21,17
%
%    Output: ns = number of states (not multiplets) in Q.
%
% Wb,Aug11,16

% outsourced from update_psi_2site.m // Wb,Dec20,21

  if nargin<2, stol=0; end
  if nargin<3, NPsi=get_NPsi_bond(Q); end

  nrm=normQS(Q);
  if ~isobject(Q), Q=QSpace(Q); end

  if ~NPsi
     if nrm>stol, Q=Q/nrm; ns=1; else Q=QSpace; ns=0; end
  elseif ~nrm, Q
     wbdie('invalid usage (got empty Q) !?');
  else

   % X=contract(Q,'!3*',Q); % overlap matrix
   % [x,I]=eigQS(X,'Rtrunc',reps); nrm=sqrt(trace(X));
   % K=QSpace(I.EK); % Q_=Q;
   % for i=1:numel(K.data)
   %     K.data{i} = diag(1./sqrt(abs(K.data{i})));
   % end
   % Q=contract(Q,QSpace(I.AK)*K);

     U=svdQS(Q,3,'stol',stol);
     Q=QSpace(U);

     if Q, d=getDimQS(Q); ns=d(end);
     else ns=0; end
  end
end

