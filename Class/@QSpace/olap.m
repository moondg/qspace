function x=olap(varargin)
% function x=olap(A [,B,opts])
%
%    calculate frobenius norm / overlap.
%
% Wb,May27,11

  getopt('init',varargin);
     if getopt('-nAB'); normA=1; normB=1; else
       normA=getopt('-nA');
       normB=getopt('-nB');
     end
  varargin=getopt('get_remaining'); narg=length(varargin);

  if narg==1, A=varargin{1}; B=A;
  elseif narg==2, A=varargin{1}; B=varargin{2};
  else varargin, error('Wb:ERR','\n   ERR invalid usage'); end
  if numel(A)~=1 || numel(B)~=1
  error('Wb:ERR','\n   ERR invalid usage'); end

  x=getscalar(QSpace(contractQS(A,1:2,B,1:2)));

  if normA, x=x/normQS(A); end
  if normB, x=x/normQS(B); end

end

