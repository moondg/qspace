function lh=scaley(sc,varargin)
% scaley - change current units for y-axis
% Usage: lh=scaley(sc, OPTS)
%
%    y -> y*sc (except if ulabel is specified: y->y/sc)
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

  getopt ('init', varargin);
     grl  = getopt('glines', []);
     ylbl = getopt('ulabel',''); if isempty(ylbl)
     tfac = getopt('-top'); else tfac=0; end
  getopt('check_error');

  h=findall(gca,'tag','scaley');
  if ~isempty(h), if numel(h)>1, wberr(...
     'got several scaley handles (%g)',numel(h)); end
     sc0=getuser(h(1),'sc'); if isempty(sc0), sc0=1; end
     delete(h);
  else sc0=1; end

  if ~isempty(ylbl)
     th=get(gca,'YLabel');
     set(th, 'String', sprintf('%s (%s)', get(th,'string'), ylbl))
     sc=1/sc;
  elseif tfac
     h=postext(-0.02,1.20,[ '\\times' num2tex(sc*sc0,'%E')],...
       'HorizontalAl','center','tag','scaley');
     setuser(h,'sc',sc); sc=1/sc;
  end

  if isequal(get(gca,'YTickMode'),'manual')
     yt=get(gca,'YTick');
     set(gca,'YTick',yt*sc);
  end
  set(gca,'YLim',get(gca,'YLim')*sc);

  for lh=[ findall(gca,'type','line'); findall(gca,'type','patch') ]'
     set(lh,'YData', get(lh,'YData')*sc);
  end

  for h=[ findall(gca,'type','text') ]'
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

