function q=using_mem(HAM)
% function q=using_mem(HAM)
%
%    Whether DMRG data is store in memory (HAM.store)
%    or mat-files (HAM.mat).
%
% Wb,Jun29,23

  q=[ isempty(HAM.store), isempty(HAM.mat) ];
  if ~xor(q(1),q(2)), wbdie('invalid storage specification'); end
  q=~q(1);

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %
