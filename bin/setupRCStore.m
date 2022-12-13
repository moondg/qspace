%  Usage #1: setupRCStore(sym [,opts])
%
%      Setup CGC container for given symmetry,
%      with the data stored in directory specified by the
%      environmental variable $RC_STORE.
%
%  Usage #2: setupRCStore(sym,q1,q2 [,opts])
%
%      explicitely generate tensor product decomposition q1*q2
%
%  Options
%
%     'npass',... number of passes through gCG data (3).
%     'dmax',...  specify max. average multiplet dimension sqrt(d1*d2)
%                 in tensor product decomposition (10*qlen).
%     '-a'        generate all, yet compatible with dmax (by default,
%                 multiplets with 2-digit qlabels are no longer used
%                 in decomposition)
%
%  (C) Wb,Aug07,14 ; Wb,Apr01,15
