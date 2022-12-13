function q=usingFullMPO(HAM)
% function q=usingFullMPO(HAM)
% Wb,Dec06,20

  if isstruct(HAM.mpo)
       q=isfield(HAM.mpo,'Q') && isfield(HAM.mpo,'data');
  else q=isa(HAM.mpo,'QSpace'); % && isfield(HAM.info,'HS')
  end

end

