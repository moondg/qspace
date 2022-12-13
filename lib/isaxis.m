function ia=isaxis(ah)
% Functions: isaxis(h)
% 
%    Check whether given input represents valid axis handles
% 
% Wb,Nov12,06

  if nargin~=1, eval(['help ' mfilename]); return; end
  if isempty(ah), ia=[]; return; end

  q=reshape(ishandle(ah),1,[]);
  for i=find(q)
     if ~isequal(get(ah(i),'Type'),'axes'), q(i)=0; end
  end

  if all(q), ia=1; 
  elseif all(~q), ia=0; else ia=reshape(q,size(ah)); end

end

