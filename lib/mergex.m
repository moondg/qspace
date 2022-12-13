function ah = mergex (varargin)
% mergex - merge axis
% Usage: ahnew = mergex (ah1,ah2,...)
%
%    extends axes ah1 to range covered by all input axes sets
%    all axes sets other than the first one in the list will
%    be deleted.
%
% Wb,Jan23,03

  ah=cat(1,varargin{:}); ah=reshape(ah,1,[]); na=length(ah);

  if ~nargin || ~all(isaxis(ah))
  eval(['help ' mfilename]); return; end

  pp=get(ah,'Position');
  pp=cat(1,pp{:});

  pp(:,3:4)=pp(:,1:2)+pp(:,3:4);

  pos=[ min(pp(:,1)), min(pp(:,2)), ...
        max(pp(:,3)), max(pp(:,4)) ];

  pos=[ pos(1:2), pos(3:4)-pos(1:2) ];

  delete(ah(2:end)); ah=ah(1);
  set(ah,'Position', pos);

end

