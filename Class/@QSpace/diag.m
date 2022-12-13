function A=diag(A,varargin)
% Function: A=diag(A [,OPTS])
%
%    extract diagonal of operator in QSpace format
%    alternatively, if only diagonal is stored,
%    diagonals become full again.
%
% Options
%
% Wb,Sep08,06

  getopt('init',varargin);
   % lflag=getopt('-l'); % lenient (accept Psi with all-in arrows) // Wb,Oct08,19
     dflag=getopt('-d');
     trans=getopt('-t');
  getopt('check_error');

  if nargin<1 || nargout>1
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

  if dflag
     if numel(A)~=1, error('Wb:ERR',...
        'Invalid QSpace (need single QSpace object)'); end
     A=diag_1(A,dflag,trans);
  else
     for k=1:numel(A)
        A(k)=diag_1(A(k),dflag,trans);
     end
  end
end

% -------------------------------------------------------------------- %
function A=diag_1(A,dflag,trans)

  if isempty(A.Q)
     if isempty(A.data)
        if dflag, A=[]; end
     else
        if numel(A.data)~=1
           error('Wb:ERR','\n   ERR unexpected QSpace structure'); end
        q=diag(A.data{1}); if trans, q=q.'; end
        if dflag, A=q; else A.data{1}=q; end
     end
     return
  end

  if numel(A.Q)~=2
     wberr('%s requires rank-2 object (%d).',mfilename,length(A.Q)); 
  elseif ~isdual_(A.Q{1},A.Q{2})
     wberr('%s requires block-diagonal operator.',mfilename);
  end

  n=length(A.data);
  if trans
       for i=1:n, A.data{i}=diag(A.data{i}).'; end
  else for i=1:n, A.data{i}=diag(A.data{i})  ; end
  end

  if dflag
     if trans, i=2; else i=1; end
     A=cat(i,A.data{:});
  end

end

% -------------------------------------------------------------------- %

