function A=mapQ(A,varargin)
% QM = Q map that transforms Q numbers to `rotated' set
% Example: SIAM: (N1,N2) --> (Q,Sz)
%
%      QM      *  Q_old   +  offset  =  Q_new
%                           
%    [ 1  1    *  [ N1    +   [-1    =  [ Q
%      1 -1 ]       N2 ]        0 ]       Sz ]
%
% Wb,Aug10,06

  if nargin<2 || nargin>3, eval(['help ' mfilename]); return; end

  for i=1:prod(size(A))
  A(i)=mapQ_1(A(i),varargin{:}); end

return

% -------------------------------------------------------------------- %

function A=mapQ_1(A,QM,offset)

  if nargin==2
     for i=1:length(A.Q)
         A.Q{i} = (QM * A.Q{i}.').'
     end
  else
     for i=1:length(A.Q)
         A.Q{i} = (QM * A.Q{i}.').' + ones(size(A.Q{i},1),1) * offset(:).';
     end
  end

return

% -------------------------------------------------------------------- %

