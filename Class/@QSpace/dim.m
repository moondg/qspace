function s=dim(A,varargin)
% function s=dim(A [,opts])
%
% Options
%
%   '-a'  all (simply return output from getDimQS())
%   '-f'  full dimension
%         (in case of CGC spaces, corresponds to state space dimension)
%   '-op' IROP dimension (i.e. checking for rank-3 QSpace)
%
% Wb,Jul31,12

  getopt('init',varargin);
     oflag=getopt('-op'); if ~oflag
     aflag=getopt('-a');
     fflag=getopt('-f'); end
  k=getopt('get_last',[]);

  if ischar(k)
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  s=getDimQS(A);

  if oflag
     r=numel(A.Q);
     if r<2 || r>3, wbdie('unexpected rank for oprator'); end
     if r<3, s=1; else; s=s(end); end
  else
     if ~isempty(k), s(:,end+1:max(k))=1; s=s(:,k); end
     if ~aflag
        if fflag, s=s(end,:); else s=s(1,:); end
     end
  end

end

