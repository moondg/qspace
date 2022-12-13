function x=subsref(G,S)
% function x=subsref(G,S)
% Wb,Nov17,11

 % x=cell(numel(G)); nG=numel(G);
 % for i=1:nG, x{i}=builtin('subsref',G(i),S); end

   x=builtin('subsref',G,S);

end

