function iz=idxzops(G)
% function iz=idxzops(G)
% Wb,Nov28,11

  iz=[]; w=0; n=numel(G);
  for i=1:n, t=G(i).type;
     if isempty(t) || numel(t)~=1, w=w+1;
     elseif isequal(t,'z'), iz(end+1)=i; end
  end

  if w, wblog('WRN',...
    'got missing or invalid operator types (%g/%g)',w,numel(G));
  end

end

