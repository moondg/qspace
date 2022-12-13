function tt=trace(G)
% Wb,Nov17,11

  tt=zeros(size(G)); nG=numel(G);
  for i=1:nG, tt(i)=trace(G(i).op); end

end

