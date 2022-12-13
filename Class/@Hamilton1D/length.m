function [L,ns]=length(HAM)
% function [L,ns]=length(HAM)
% Wb,Sep29,15

  L=numel(HAM.mpo);
  if nargout>1, ns=unique([HAM.mpo.stype]); end

end

