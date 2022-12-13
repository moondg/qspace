function i=finditag(A,pat,nx)
% function i=finditag(A,pat [,n])
%
%    find itag with specit regexp pattern
%    if n is specified, the result must yield exactly n matches.
%
% Wb,Jun03,19

  if nargin<2 || numel(A)~=1
     helpthis, if nargin || nargout, wberr('invalid usage'), end
     return
  end

  if isempty(A.info), i=[];
  else i=find(cellfun(@isempty,regexp(A.info.itags,pat))==0);
  end

  if nargin>2
     if nx<1, wberr('invalid usage (enforcing n=%d !?)',nx); end
     n=numel(i); if n~=nx
        if ~n, wberr('failed to find match for ''%s''',pat);
        elseif nx==1
             wberr('got multiple matches for ''%s'' (%d)',pat,n);
        else wberr('unexpected number of matches for ''%s'' (%d/%d)',pat,n,nx);
        end
     end
  end

end

