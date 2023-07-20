function xx=getxlim(varargin)
% Function xx=getxlim([OPTS])
% Options
%
%    '-view'   determine XLim fully shows all y-data in given view
%    '-data'   default: xlim is full x-data range
%
% Wb,Jan10,08 - see also xtight.m Wb,2002.

% if nargin>1
%    eval(['help ' mfilename]);
%    if nargin || nargout, wbdie('invalid usage'), end, return
% end

  o={}; if isempty(varargin), varargin{1}='-data'; end

  switch varargin{1}
    case '-pos',  o={'-pos'};
    case '-view', xx=getXLimView;
    case '-data', xx=getXLimData(o{:});
    otherwise
       eval(['help ' mfilename]);
       wbdie('invalid usage')
  end

end

% -------------------------------------------------------------------- %
function xx=getXLimView()

  lh=findall(gca,'Type','Line','visible','on');
  m=length(lh); xx=zeros(m,3);

  yl=ylim;

  islog=isequal(get(gca,'XScale'),'log');

  for i=1:m
      yd=get(lh(i),'YData'); ii=find(yd>=yl(1) & yd<=yl(2));
      xd=get(lh(i),'XData'); xd=xd(ii);
      if islog, xd=xd(find(xd>0)); end
      if ~isempty(xd)
           xx(i,:)=[min(xd), max(xd), length(yd)];
      else xx(i,:)=[+1 -1  0]*1E99;
      end
  end

  i=find(xx(:,3)>2); if isempty(i), i=1:size(xx,1); end
  xx=[ min(xx(i,1)), max(xx(i,2)) ];

  if xx(1)==1E99
     wblog('WRN','failed to determine XLim for data within view');
     xx=xlim;
  end

end

% -------------------------------------------------------------------- %
function xx=getXLimData(varargin)

  if nargin
     if ~isequal(varargin{1},'-pos'), wbdie('invalid usage'); end
  pflag=1; else pflag=0; end

  lh=findall(gca,'Type','Line','visible','on');
  m=length(lh);

  mark=zeros(size(lh));
  for i=1:m
     if ~isempty(regexp(get(lh(i),'tag'),'mark')), mark(i)=1; end
  end
  lh(find(mark))=[]; m=length(lh); xx=zeros(m,2);

  for i=1:m
      xd=get(lh(i),'XData'); if pflag, xd=xd(find(xd>0)); end
      if ~isempty(xd)
           xx(i,:)=[min(xd), max(xd)];
      else xx(i,:)=[+inf -inf];
      end
  end

  xx=[ min(xx(:,1)), max(xx(:,2)) ];

end

% -------------------------------------------------------------------- %

