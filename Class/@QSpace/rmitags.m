function A=rmitags(A)
% function A=rmitags(A)
%
%    skip itags of given QSpace(s)
%
% *** NB! deprecated; rather use untag ***
% Wb,Apr14,16

% tags: skipitags, skip_itags

  if nargin>1
     helpthis, if nargin || nargout
     wbdie('invalid usage'), end, return
  end

  for k=1:numel(A)
     A(k).info.itags=regexprep(A(k).info.itags,'[^\*]*','');
  end

end

