function hh=splitax(ah,varargin)
% splitax - split axis
%
% Usage: ahs = splitax(ah, [, [n m], parameters])
%        splits axes ah into nxm (default: 1x2)
%        upper left axes will contain orginial ah axes set
%
% Parameters:
%
%    'tight' no margin in between
%    'dx'    margin in between the two new axes (in points; 50)
%    'dy'    margin in between the two new axes (in points; 50)
%
% Remaining parameters are applied to all axes handles (old+new).
%
% Wb,Feb05,03  Wb,Mar14,07

  if nargin<1 || ~isscalar(ah) || ~isaxis(ah)
     helpthis, if nargin || nargout, wberr('invalid usage'), end
     return
  end

  if nargin>1 && isnumeric(varargin{1}) && isequal(size(varargin{1}),[1 2])
     m=varargin{1}(1);
     n=varargin{1}(2); varargin(1)=[];
  else
     m=1; n=2;
  end

  dx=50; dy=50;

  getopt ('init', varargin);
     tight = getopt('tight'); if tight, dx=0; dy=0; end
     dx_in = getopt('dx',[]);
     dy_in = getopt('dy',[]);
  varargin=getopt('get_remaining');

  U0=get(ah,'Units');
  set(ah,'Units','points');
  p0=get(ah,'Position');

  if ~isempty(dx_in), dx=dx_in;
  elseif dx>0
     if (p0(3)+dx)/n-dx<2*dx
     dx=(p0(3)+dx)/n * 0.16; end
  end

  if ~isempty(dy_in), dy=dy_in;
  elseif dy>0
     if (p0(4)+dy)/m-dy<2*dy
     dy=(p0(4)+dy)/m * 0.20; end
  end

  rs(1)=(p0(3)+dx)/n-dx;
  rs(2)=(p0(4)+dy)/m-dy;

  if rs(1)<dx || rs(2)<dy, wberr([
    'axis too small to be split any further for requested margin' ...
    'dxy=(%g,%g) having (%g,%g)\n'],dx,dy,rs(1),rs(2));
  end

  p0=[p0(1:2), rs]; hh=zeros(m,n);

  rs=rs+[dx dy];

  for i=1:m
  for j=1:n
     pp=p0+[(j-1)*rs(1), (m-i)*rs(2), 0, 0];

     if i>1 || j>1
          hh(i,j)=axes('Units','points','Position', pp);
     else hh(i,j)=ah; set(hh(1),'Units','points','Position', pp);
     end
  end
  end

  set(hh(2:end),'box',get(hh(1),'box'),...
    'LineW', get(hh(1),'LineW'), ...
    'FontSize', get(hh(1),'FontSize'));

  set(hh,'Units',U0);

  if ~isempty(varargin), set(hh,varargin{:}); end

  if tight
     set(hh(1:end-1,:),'XTickLabel', []);
     set(hh(:,2:end),  'YTickLabel', []);
  end

end

