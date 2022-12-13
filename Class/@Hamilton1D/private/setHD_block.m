function H=setHD_block(H,k,l,h)
% function H=setHD_block(H,k,l,h)
% set (block) matrix elements H(k,l)=h [and H(l,k)=h' if k~=l]
% Wb,Mar24,16

% outsourced from update_psi_2site.m // Wb,Dec20,21

  if numel(H.Q)~=2, H
     wberr('unexpected H struct !?'); 
  end
  if isempty(h.data)
     return
  end

  if isempty(h.Q), h=getscalar(h);
     H.data{1}{k,l}=h;  if k~=l
     H.data{1}{l,k}=h'; end
     return
  end

  if numel(H.Q)~=2 || numel(h.Q)~=2, H, h
     wberr('unexpected objects !?'); 
  end

  [i1,i2,I]=matchIndex(H.Q{1},h.Q{1});
  if ~isempty(I.ix2) || numel(i2)~=size(h.Q{1},1), H, h
     wberr('got new symmetry sectors in h !?'); end

  xflag=(min(k,l)==1);

  for i=1:numel(i2), x=h.data{i2(i)};
      H.data{i1(i)}{k,l}=x;  if k~=l
      H.data{i1(i)}{l,k}=x'; end
  end

end

