function r=symrank(A)
% function r=symrank(A)
% Wb,Jan09,17

% see also $SYM/private/get_rank.m

  if numel(A)~=1
     helpthis, error('Wb:ERR','invalid usage')
  end

  r=getsym(A,'-d');

end

