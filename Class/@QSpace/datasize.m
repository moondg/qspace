function dd=datasize(A,sflag)
% Function [D,DD]=datasize(A [,'-s'])
%
%    get full listing of size(data{}),
%    except if option '-s' is specified which gives summary 
%    in terms of total size in bytes of A.data
%
% Wb,Apr24,08

  if ~nargin
     helpthis, if nargin || nargout, wberr('invalid usage'), end
     return
  end

  if nargin<2, n=length(A.data); r=length(A.Q); dd=ones(n,r);
     for i=1:n
        s=size(A.data{i});
        dd(i,1:length(s))=s;
     end
  else
     if ~isequal(sflag,'-s'), wberr('invalid usage (2nd arg)'); end
     for i=numel(A):-1:1, q{i}=A(i).data; end
     s=whos('q'); dd=s.bytes;
  end

end

