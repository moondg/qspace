function [r,nq]=symrank(A)
% function [r,nq]=symrank(A)
%
%    Get symmetry rank for each symmetry.
%
% Wb,Jan09,17

% see also $SYM/private/get_rank.m

  if numel(A)~=1
     helpthis, wbdie('invalid usage')
  end

  r=getsym(A,'-d');
  if nargout>1
     if isempty(A.Q), nq=0; else nq=size(A.Q{1},2); end
     if nq~=sum(r)
        wbdie('mismatch in number of symmetry labels (%d/%d)',nq,sum(r));
     end
  end

end

