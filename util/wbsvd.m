%  Usage: [U,S,V] = wbsvd(H [,I]);
%
%      singular value decomposition based on BLAS dgesvd
%      for genearlized rank-r object.
%      index set I specifies which dimensions are combined
%      into first fused dimension. For rank-2 objects
%      and I omitted, this corresponds to the regular SVD
%      on matrizes.
%
%  AW (C) Sep 2009
