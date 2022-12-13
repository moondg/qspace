function x=olap(A,B,varargin)
% function x=olap(A,B)
%   
%    calculate overlap between A and B as in conj(A(:))'*B(:).
%    for vectors this is usual <A|B>, while for matrizes
%    this corresponds to trace(A'*B).
%
% Options
%
%    '-nA'  normalize with respect to A (e.g. devide by <A|A>)
%    '-nB'  same as -nA but using B, instead.
%
% Wb,Jun08,10

  getopt('init',varargin);
     normA=getopt('-nA'); if ~normA,
     normB=getopt('-nB'); else normB=0; end
  getopt('check_error');

  x=conj(A(:)'*B(:));

  if normA, q=norm(A,'fro')^-2; x=x*q; end
  if normB, q=norm(B,'fro')^-2; x=x*q; end

end
