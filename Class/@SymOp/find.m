function i=find(CC)
% function i=find(CC)
% Wb,Sep27,12

  mark=ones(size(CC)); n=numel(CC); s={};
  for i=1:n
     if isempty(CC(i).op), mark(i)=0; continue; end
     if ~isempty(s)
        if ~isequal(size(CC(i).op),s)
        error('Wb:ERR','\n   ERR incompatible SymOp set'); end
     else
        s=size(CC(i).op);
     end
  end

  i=find(mark);

end

