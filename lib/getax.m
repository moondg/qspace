function ah=getax(varargin)
% Function getax(atag)
%
%    Get axes set with given tag.
%    See also setca.m
%
% Wb,Dec03,07

  if length(varargin)~=1 || ~ischar(varargin{1})
     eval(['help ' mfilename]);
     if nargin || nargout, wbdie('invalid usage'), end, return
  end

  atag=varargin{1};
  ah=findall(groot,'type','axes','tag',atag);

  if isempty(ah)
     wblog(2,'No axes with tag `%s'' found.',atag);
     return
  end

  if length(ah)>1
     h=findall(gcf,'type','axes','tag',atag);
     if length(h)==1, ah=h; else
        wblog(2,'multiple axes sets with tag `%s'' found (%d)',atag,length(ah));
        ah=ah(1);
     end
  end

end

