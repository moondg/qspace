function q=using_full_MPO(HAM)
% function q=using_full_MPO(HAM)
% Wb,Dec06,20

  q=false;
  if isfield(HAM,'mpo')
     if isstruct(HAM.mpo)
          q=isfield(HAM.mpo,'Q') && isfield(HAM.mpo,'data');
     else q=isa(HAM.mpo,'QSpace');
     end
  end

end

