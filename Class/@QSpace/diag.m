function A=diag(A,varargin)
% Function: A=diag(A [,OPTS])
%
%    Toggle that switches between diagonal formats similar to matlab's
%    behavior of diag(): for full operator A, extract diagonal of A
%    in QSpace format; if already in diagonal format, diagonals become
%    restored to full format again.
%
% Options
%
%    -d    return numeric data vector of diagonal only
%    -t    store diagonal as rows (by default,
%          diag(M) of matrix M puts diagonal of M into column)
%    -c    compact diagonal QSpace to diagonal only
%          this throws an error if non-diagonal data is present in A
%
% Wb,Sep08,06

  getopt('init',varargin);
   % lflag=getopt('-l'); % lenient (accept Psi with all-in arrows) // Wb,Oct08,19
     dflag=getopt('-d');
     trans=getopt('-t');
     cflag=getopt('-c');
  getopt('check_error');

  if nargin<1 || nargout>1
     eval(['help ' mfilename]);
     if nargin || nargout, wbdie('invalid usage'), end, return
  end

  if dflag
     if numel(A)~=1, wbdie(...
       'invalid usage (single QSpace object required)'); end
     A=diag_1(A,cflag,dflag,trans);
  else
     for k=1:numel(A)
        A(k)=diag_1(A(k),cflag,dflag,trans);
     end
  end
end

% -------------------------------------------------------------------- %
function A=diag_1(A,cflag,dflag,trans)

  nd=numel(A.data);

  if cflag || dflag, fullD=zeros(1,nd);
     for i=1:nd
        s=size(A.data{i}); n=prod(s);
           if numel(s)~=2, wbdie('unexpected data size'); end
           if n==1, continue; end
        j=find(s~=1); m=numel(j);

        if m==1, fullD(i)=-j;
        else fullD(i)=2;
           if cflag, M=A.data{i};
              j=1:min(s); j = j + s(1)*(j-1);
              x=norm(M(j)); M(j)=0;
              e=norm(M,'fro')/max(1,x); if e>1E-12
              wbdie('got non-diagonal data @ %.3g (got option -c)',e); end
           end
        end
     end
     q=unique(fullD(find(fullD)));
     if isempty(q), q=1;
     elseif numel(q)>1
        wbdie('invalid usage (got mixed diagonal setting)'); 
     end
     fullD=q;
  end

  if isempty(A.Q)
     if ~nd
        if dflag, A=[]; end
     else
        if nd~=1, wbdie('unexpected QSpace structure'); end
        q=diag(A.data{1}); if trans, q=q.'; end
        if dflag, A=q; else A.data{1}=q; end
     end
     return
  end

  if numel(A.Q)~=2
     wbdie('%s requires rank-2 object (%d)',mfilename,length(A.Q)); 
  elseif ~isdual_(A.Q{1},A.Q{2})
     wbdie('%s requires block-diagonal operator',mfilename);
  end

  if dflag, A=A.data;
       if fullD>0, for i=1:nd, A{i}=diag(A{i}); end
            A=cat(1,A{:});
       else A=cat(-fullD,A{:}); end
       if trans, A=A.'; end
  elseif trans
       for i=1:nd, A.data{i}=diag(A.data{i}).'; end
  else for i=1:nd, A.data{i}=diag(A.data{i})  ; end
  end

end

% -------------------------------------------------------------------- %

