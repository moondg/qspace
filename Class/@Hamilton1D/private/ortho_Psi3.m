function [Q,ns,nrm,msg]=ortho_Psi3(Q,stol,NPsi,kflag)
% function [Q,ns,nrm,msg]=ortho_Psi3(Q,stol,NPsi, [kflag])
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

  if nargin<4, kflag=0;
  elseif isequal(kflag,'-k'), kflag=1; 
  elseif isequal(kflag,'-K'), kflag=2; 
  else kflag, wblog('WRN','got unexpected kflag (ignore)'); kflag=0;
  end

  nrm=normQS(Q); msg='';
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

     if kflag<=0, Q=QSpace(U);
     else 
      % NB! SVD Q -> U does not stay close to the original basis
      % => orthonormal basis needs to stay close // Wb,Apr16,24
      %    e.g. for correpondence Eg <> Psi index in Rho
      %    see MAT/matrices.tex -> Closest orthonormal basis
      % 1) project to arbitrary but fixed orthonormal basis U*U' (see also 3)
        X=contractQS(U,'!3*',Q); Q_=Q;
        for i=1:numel(X.data)
           [u,s,v]=svd(X.data{i});
           X.data{i}=u*v';
        end
        Q=QSpace(contractQS(U,3,X,1));

        e=norm(Q-Q_);
        if e>5e-4
           msg=sprintf('δψ=%.1e',e);
           if e>0.05
              if kflag>1 || nargout<4
                   wblog('WRN', msg); if kflag>1, wbstop, end
              else msg=[ 'WRN ' msg]; end
           elseif nargout<4
              wblog(' * ', msg);
           end
        end
     end

     if Q, d=getDimQS(Q); ns=d(end);
     else ns=0; end
  end
end

