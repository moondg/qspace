function ia=isfig(hh)
% Functions: isfig(h)
% 
%    Check whether given input is valid figure handles
% 
% Wb,Nov12,06

  if nargin~=1, eval(['help ' mfilename]); return; end

  if ~ishandle(hh), ia=0; return; end

  m=numel(hh); ia=zeros(size(hh));
  for i=1:m
     if ishandle(hh(i)) && isequal(get(hh(i),'Type'),'figure')
     ia(i)=1; end
  end

end

