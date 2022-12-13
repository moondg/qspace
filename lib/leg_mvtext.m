function hh=leg_mvtext(l0,dy)
% Function: leg_mvtext(l0,dy)
%
%    with dy either a scalar (for all legend entries
%    or a vector for every legend entry
%
% Wb,May07,07

  if nargin~=2 || numel(l0)~=1 || ~ishandle(l0)
     helpthis, if nargin || nargout, wberr('invalid usage'), end
     return
  end

  t=get(l0,'type');
  if isequal(t,'axes')
     hh=findall(l0,'type','text','tag','')';
  elseif isequal(t,'legend')
     hh=get(l0,'ItemText');
  else l0, t, wberr('invalid handle'); end

  m=numel(hh);

  if length(dy)==1, dy=dy(ones(m,1));
  end
  n=min(numel(hh),numel(dy));

  for i=1:n, h=hh(end-i+1);
     set(h,'Pos',get(h,'Pos')+[0 dy(i) 0]);
     setuser(h,'leg_dy',dy(i));
  end

  if ~nargout, clear hh; end

end

