function i=sameas(A,B,varargin)
% function i=sameas(A,B [,eps=1E-10][,'-l'])
%
%    check whether QSpaces A and B are the same
%    except for minor differences such as 
%    info.cgr.cid(end) and info.cgr.nnz
%
% Options
%
%    eps    data is the same to within eps (default: 1E-10)
%    -l     permit different order entry (default: no)
%
% Wb,Aug08,16

  i=1;
  if ~isequal(size(A),size(B)), i=0; return; end

  eps=1E-10; lflag=0;
  if nargin==3
     if isnumeric(varargin{1}) && abs(varargin{1})<1, eps=varargin{1};
     else lflag=varargin{1}; end
  elseif nargin==4, [eps,lflag]=deal(varargin);
  elseif nargin>4, wbdie('invalid usage'); end

  if ~isequal(lflag,0)
     if isequal(lflag,'-l'), lflag=1; 
     elseif ~isnumeric(lflag), wbdie('invalid usage'); end
  end

  if ~isa(B,'QSpace'), B=QSpace(B); end

  for k=1:numel(A)
     if ~isequal(A(k),B(k)), Ak=A(k);
        if ~isequal(Ak.Q,B(k).Q), if ~lflag, i=0; break; end
           if ~isequal(size(Ak.Q),size(B(k).Q)) || ~isempty(Ak.Q) && ...
              ~isequal(size(Ak.Q{1}),size(B(k).Q{1})); i=0; break; end

           QA=A(k).Q; QA=[QA{:}];
           QB=B(k).Q; QB=[QB{:}]; [ib,ia,Im]=matchIndex(QB,QA,'-s');
           if ~isempty(Im.ix1) || ~isempty(Im.ix2), i=0; break; end
           if ~isequal(ib,1:size(QB,1)), i=0; break; end
           Ak=getsub(A(k),ia);
           if ~isequal(Ak.Q,B(k).Q), i=0; break; end
        end

        if ~isequal(getqdir(Ak), getqdir(B(k))), i=0; break; end

        e=[normQS(Ak), normQS(B(k))];
        if any(e>eps)
           e=normQS(Ak-B(k))/max(e);
           if e>eps, i=0; break; end
        end
     end
  end

end

