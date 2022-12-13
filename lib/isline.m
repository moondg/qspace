function r=isline(hh)
% Functions: i=isline(h)
% 
%    Check whether given input is valid line handles
% 
% Wb,Nov12,06

  r=0;
  if ~nargin, eval(['help ' mfilename]); return; end
  hh=reshape(hh,[],1);

  if any(~ishandle(hh)), return; end

  tt=get(hh,'Type'); if ~iscell(tt), tt={tt}; end
  if any(cellfun(@isempty, ...
     regexpi(tt,'line|errorbar|hggroup')))
     return
  end

  r=1;
end

