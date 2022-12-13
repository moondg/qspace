function G=untrace(G)
% function G=untrace(G)
% make operators G traceless
% Wb,Nov18,11

  nG=numel(G);
  for i=1:nG, q=G(i).op; t=trace(q);
     if abs(t)>1E-12, n=size(q,1);
        G(i).op=q-diag(repmat(t/n,1,n));
     end
  end

end

