function q=qmat2char(q)
% function qs=qmat2char(q)
% See also char2qmat.m
% Wb,Aug18,16

% see also $MLIB/int2str2.m -> int2str2(q,'-X') // Wb,Oct02,19

  if any(reshape(q~=round(q),[],1)), q
     wbdie('invalid usage (assuming integer q-labels) !?'); 
  end

  if any(reshape(q<0,[],1))
     if any(reshape(q>9,[],1)) || size(q,1)~=1, q
        wbdie('invalid usage (assuming non-negative ints) !?'); 
     end
     i=find(q<0); q=char(abs(q)+'0');
     q(end+1,:)=' ';
     q(end,i)='-';
  else
     i=find(q>9); q(i)=q(i)+7; q=char(q+'0');
  end

end
