function x=mpsoverlap(A,varargin)
% function x=mpsoverlap(A [,B])
%
%    calculate MPS overlap <B|A> (note the reversed order!)
%    if B is not specified, B:=A, by default.
%
% adapted from VER/vbsoverlap.m

  getopt('init',varargin);
     qflag=getopt('-q');
  B=getopt('get_last',A);

  x=1; N=numel(A);

  s=[ size(A{1},1), size(B{1},1) ];
  if any(s~=1)
     if ~norm(diff(s)), x=eye(s); if ~qflag,
        wblog('WRN','L-vac has dim=%g - assume L-orthogonal',max(s)); end
     else
        error('Wb:ERR','incompatible L-vac in A and B (dim=%g,%g)',s);
     end
  end

  for k=1:N,
      x=contract(contractmat(A{k},x,1), conj(B{k}),[1 3],[1 3]);
  end

  s=size(x);
  if any(s~=1)
     if ~qflag
        wblog('WRN','R-vac has dim=%g - assume R-orthogonal',max(s)); end
     x=trace(x);
  end

end

