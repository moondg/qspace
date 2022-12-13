function mv2front(lh)
% Function - mv2front
%
%    move given (line) handles of current axis to foreground
%
% Usage: mv2front(lh)
%
% keywords: get line
% Wb,Jun08,06

  ah=get(gca,'children');

  ia=[]; ie=[];
  for i=1:length(lh)
     k=find(ah==lh(i));
     if ~isempty(k), ia=[ia,k]; else ie=[ie,i]; end
  end

  if ~isempty(ie)
     wblog('ERR Some handles missing in current axes (%d)!', length(ie));
     lh(ie)=[];
  end

  ah(ia)=[];

  set(gca, 'children', [lh(:); ah(:)]);

return

