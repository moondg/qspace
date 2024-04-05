function lh=scaley(varargin)
% function lh=scaley([ah,] sc [,opts])
%
%    change current units of y-axis y -> y*sc 
%    (except if ulabel is specified, then y->y/sc)
%    for axis ah (default: gca)
%
% Options:
%
%    'glines', ...  draw vertical grid lines at positions indicated,
%                   with grl in units of sc
%    'ulabel',...   unit label to use, as in <current ylabel> (<ulabel>)
%    '-top'         without any further argument puts sc factor at the
%                   left top of the axis set (as in \times 10^{-3}) if sc=1E-3
%
% Wb,Oct24,05

  if nargin && isaxis(varargin{1}), l=2;
     ah=varargin{1}; n=numel(ah);
     if n~=1, wbdie('invalid usage (got %d axis handles)',n); end
  else ah=gca; l=1;
  end

  if nargin>=l, sc=varargin{2}; l=l+1;
  else wbdie('invalid usage'); end

  getopt('init',varargin(l:end));
     grl  = getopt('glines', []);
     ylbl = getopt('ulabel',''); if isempty(ylbl)
     tfac = getopt('-top'); else tfac=0; end
  getopt('check_error');

  h=findall(ah,'tag','scaley');
  if ~isempty(h), if numel(h)>1, wbdie(...
     'got several scaley handles (%g)',numel(h)); end
     sc0=getuser(h(1),'sc'); if isempty(sc0), sc0=1; end
     delete(h);
  else sc0=1; end

  if ~isempty(ylbl)
     th=get(ah,'YLabel');
     set(th, 'String', sprintf('%s (%s)', get(th,'string'), ylbl))
     sc=1/sc;
  elseif tfac
     h=postext(-0.02,1.20,[ '\\times' num2tex(sc*sc0,'%E')],...
       'HorizontalAl','center','tag','scaley');
     setuser(h,'sc',sc); sc=1/sc;
  end

  if isequal(get(ah,'YTickMode'),'manual')
     yt=get(ah,'YTick');
     set(ah,'YTick',yt*sc);
  end
  set(ah,'YLim',get(ah,'YLim')*sc);

  for lh=[ findall(ah,'type','line'); findall(ah,'type','patch') ]'
     set(lh,'YData', get(lh,'YData')*sc);
  end

  for h=[ findall(ah,'type','text') ]'
     if ~isequal(get(h,'units'),'data') || islabel(h), continue, end
     p=get(h,'Position'); p(2)=p(2)*sc;
     set(h,'Position',p);
  end

  if ~isempty(grl)
     hold on
     n=length(grl); lh=zeros(length(n));
     for i=1:n
         lh(i)=plot([grl(i), grl(i)], ylim, 'k--');
     end
     yl=[min(grl), max(grl)]; yl = yl+0.2*[-1 1].*diff(yl); xlim(yl);
  end

  if ~nargout, clear lh; end

end

