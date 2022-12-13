function P=update_size(P)
% function P=update_size(P)
% Wb,May12,09

   P.r=numel(P.vars);

   P.s=ones(1,P.r);
   for i=1:P.r, P.s(i)=numel(P.data{i}); end

   P.n=prod(P.s);

end

