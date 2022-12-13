function q=getenvb(v,bb)
% function getenvb(var [,bb])
% Wb,Dec19,21

  q=getenv(v);
  if isempty(q), q=0; return; end
  if ~regexp(q,'^\d+$'), wblog('WRN',['unexpected ENV %s ' ...
    '(expeting integer interpreted as bit set)'],v);
     q=0; return
  end

  q=str2num(q);
  if nargin<2 || isequal(bb,'*') || isempty(bb), return; end

  if ~isnumeric(bb), wbdie('invalid usage (got non-numeric bb)');
  elseif any(bb<1)
     wbwrn('invalid usage (bit index is 1-based)'); bb(bb<1)=[];
  elseif any(bb>32), wbwrn('invalid usage (bb is bit index, not value)'); end

  q=bitand(q, sum(2.^unique(bb-1)));

end

