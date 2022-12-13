function ns=numsym(A)
% function ns=numsym(A)
%
%    determine number of symmetries in QSpace(s) A.
%
% Wb,Feb06,14

  if nargin>1
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  ns=zeros(size(A));
  for k=1:numel(A)
     if isempty(A(k).info) || isempty(A(k).info.cgr)
        if isempty(A(k).Q) || isempty(A(k).Q{1}), error('Wb:ERR',['\n   ' ... 
          'ERR failed to determine number of symmetries for empty QSpace']);
        end
        ns(k)=size(A(k).Q{1},2);
     else
        if ~isfield(A(k).info,'qtype') || isempty(A(k).info.qtype)
           error('Wb:ERR',['\n   ' ... 
          'ERR failed to determine number of symmetries for empty QSpace']);
        end
        q=strread(A(k).info.qtype,'%s','delimiter',',');
        ns(k)=numel(q);
     end
  end

end

