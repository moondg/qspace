function x=olap(A,B,varargin)
% function x=olap(A,B)
% Wb,Dec05,11

  if numel(A)~=1 || numel(B)~=1
     wbdie('invalid usage'); end

  x=olap(A.op,B.op,varargin{:});

end

