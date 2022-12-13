function is=got(P,field)
% function i=got(P,field)
% check wheter given PSet contains field P
% Wb,Jul01,11

  if ~isa(P,'PSet'), P=PSet(P); end
  nn=P.vars; n=numel(nn); is=0;
  for i=1:n
     if isequal(nn{i},field), is=1; break; end
  end

end

