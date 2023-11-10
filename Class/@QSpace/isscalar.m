function q=isscalar(A,varargin)
% function i=isscalar(A [,opts])
%
%    Checks whether QSpace is a scalar, meaning rank-0
%    with a single number in data [e.g., a scalar is obtained
%    by tracing/contracting all indices resulting in a single
%    number].
%
% Options
%
%     -l   lenient (also consider a QSpace with non-empty Q
%          a scalar, as long as data contains a single number)
%     -x   extended mode (in case of input QSpace array,
%          show isscalar() for each individual QSpace
%
% Wb,Aug04,08

  getopt('init',varargin);
     lflag=getopt('-l');
     xflag=getopt('-x');
  getopt('check_error');

  nA=numel(A); q=false(1,nA);
  for i=1:numel(A), q(i)=isscalar_1(A(i),lflag); end

  if xflag, q=reshape(q,size(A));
  else
     q=all(q);
  end

end

% -------------------------------------------------------------------- %

function q=isscalar_1(A,lflag)

  nq=numel(A.Q); if nq, nq=size(A.Q{1},1); end
  nd=numel(A.data);
  if nq && nq~=nd
     wbdie('severe QSpace inconsistency (%g/%g) !?',nq,nd);
  end

  if nd==1, nd=numel(A.data{1}); end
  if nd==1, q=true; else q=false; end
  if ~q || (q && ~nq), return; end

  if ~lflag, q=false; end

end

% -------------------------------------------------------------------- %

