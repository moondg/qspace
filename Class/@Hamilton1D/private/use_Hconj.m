function q=use_Hconj(HAM)
% function q=use_Hconj(HAM)
%
%    The full MPO may still require "H.c." as indicated by
%    info.mpo.use_hconj, i.e., in which case Hfull = H_mpo + H_mpo^\dagger
%
% Wb,Dec04,20

  if usingFullMPO(HAM)
       q=HAM.info.mpo.use_hconj;
  else q=0; end

end

