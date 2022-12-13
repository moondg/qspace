function q=qmat2char(q)
% function qs=qmat2char(q)
% See also char2qmat.m
% Wb,Aug18,16

% see also $MLIB/int2str2.m -> int2str2(q,'-X') // Wb,Oct02,19

  if any(reshape(q<0,[],1)) || any(reshape(q~=round(q),[],1)), q
     error('Wb:ERR','\n   ERR invalid usage (assuming non-negative ints) !?'); 
  end

% diff(double('9A'))=8
  i=find(q>9); q(i)=q(i)+7; q=char(q+'0');

end
