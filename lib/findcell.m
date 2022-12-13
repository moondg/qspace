function ii=findcell(c,a)
% Function: ii=findcell(c,a)
%
%    Searches cell c for entries a.
%    where c can contain a mix of objects, such as string, numbers, etc.
%
% Wb,Oct23,06 ; Wb,Jun07,19

  ii=find( cellfun(@(x) isequal(x,a), c));

end

