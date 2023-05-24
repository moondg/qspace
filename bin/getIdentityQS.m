%  Usage 1: A=getIdentityQS(A [,i1, perm,'-0']);
%
%     get plain identity operator from given QSpace.
%
%   '-0'  this generates the 1J symbol, i.e. an identity operator
%         with all indizes inward. This can be used to "revert
%         arrows", and subsequently also to avoid dragging along
%         a scalar singleton index if the total symmetry of a
%         given wave function has q=0.
%
%  Usage 2: A=getIdentityQS(A [,i1 [,at]], B [,i2 [,bt] itag, perm]);
%
%     get tensor product space of input spaces defined by
%     A and B. If i[12] is not specified, rank-2 objects
%     are assumed taking (the first) two indizes.
%
%     Default index order of the output space is (1,2,12*)
%     with the itags of A and B inherited if present.
%     An <itag> name for the combined output space (dim3) may be
%     specified (conj flag if present, will be ignored).
%
%     [06/01/2019] In addition, for an input A-tensor with LRs
%     index order convention, its itags may be inherited to C
%     by specifying <at> xor <bt> (see usage above) in the format
%     '-A[:..]' where `..' indicates extra characters to be
%     catenated at the end of the itag for the fused index.
%
%     For example, consider an A-tensor with itags {K01,K02*,s02},
%     and E a local identity operator without itags; then in order
%     to ensure a complete local state space, one may write
%     >> getIdentityQS(A,1,'-A:~',E)
%     which generates a QSpace with itags {K01,s02,K02~*}.
%     K01 is inherited from A via `A,1'; the remaining two itags
%     are derived from the remaining two itags in A;
%     therefore `A,i' must have i=1 or 2 with in LRs index order.
%     conj-flags are properly adjusted as needed.
%     Conversely, getIdentityQS(A,2,'-A:~',E)
%     generates itags {K02,s02,K01~*}, e.g. used for R->L sweep.
%
%  Finally, the output index order may be changed by specifying
%  a permutation perm to be applied on the final object.
%  If itag is not specified, perm must be in numerical format,
%  otherwise a (compact) string format is also accepted.
%
%  Further options
%
%    '-h'  display this usage and exit
%    '-v'  verbose
%
%  Note that mex files do not allow to return class objects,
%  hence A is returned as a QSpace structure. To get a class
%  object, use QSpace(getIdentityQS(...)), or equivalently,
%  if getIdentity() is properly defined as a wrapper routine
%  within matlab's Class/@QSpace (see MPS Pack), having a QSpace
%  input A, this may be shortended to getIdentity(A,...).
%
%  (C) AW : Apr 2010 ; Oct 2014 ; May 2017
