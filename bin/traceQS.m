%  Usage: X = traceQS(A [,i1,i2])
%
%      Trace out index set i1 with i2 of given QSpace object.
%      For example, with i1='12' or [1 2] and i2='34' or [3 4]
%      (also permitting compact string-based index convention
%      as with contractQS), this contracts index 1 with 3,
%      and index 2 with 4.
%
%      For tensors A of even rank r when specifying no index sets
%      this is interpreted as the full trace of indices i1=1:r/2
%      with i2=r/2+1:r.
%
%      If not all indices of A are traced out, the returned X
%      is a QSpace of rank r-2*length(i1). If all indices are
%      traced out, the returned X is a plain number.
%
%  AW (C) 2006 ; 2018 ; 2025
