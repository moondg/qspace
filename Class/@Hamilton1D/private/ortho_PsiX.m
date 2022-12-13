function [Q,X]=ortho_PsiX(Q,U,NPsi)
% function [Q,x]=ortho_PsiX(Q,U)
%
%    project Q orthogonal to U assuming lr[g] index order
%    (NB! this assumes that U does already represent an
%    orthonormal state space!)
%
% Wb,Aug11,16

% outsourced from update_psi_2site.m // Wb,Dec20,21

  X=contract(U,'12*',Q,'12');
  if ~NPsi
       Q=Q-getscalar(X)*U;
  else Q=Q-contract(U,'3',X,'1');
  end

end

