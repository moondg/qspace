function ah=setca(varargin)
% Function setca(tag)
%
%    Set axes set with given tag the current one.
%    (also sets CurrentFigure without redrawing it)
%
% Wb,Dec03,07

  getopt('init',varargin);
     dflag =getopt('-disp');
  varargin=getopt('get_remaining');

  if length(varargin)~=1 || ~ischar(varargin{1})
     eval(['help ' mfilename]);
     if nargin || nargout, wbdie('invalid usage'), end, return
  end

  tag=varargin{1};
  ah=findall(groot,'type','axes','tag',tag);

  if isempty(ah)
     wblog(2,'No axes with tag `%s'' found.',tag);
     return
  end

  if length(ah)>1
     h=findall(gcf,'type','axes','tag',tag);
     if length(h)==1, ah=h; else
        wblog(2,'multiple axes sets with tag `%s'' found (%d)',tag,length(ah));
        ah=ah(1);
     end
  end

  if ~dflag
     fh=get(ah,'Parent');
     set(0,'CurrentFigure',fh);
     set(fh,'CurrentAxes', ah);
  else
     axes(ah);
  end

  if ~nargout, clear ah; end

end

