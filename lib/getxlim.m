function xx=getxlim(varargin)
% Function xx=getxlim([ah,][opts])
% Options
%
%    '-view'   determine XLim fully shows all y-data in given view
%    '-data'   default: xlim is full x-data range
%
% Wb,2002; Wb,Jan10,08 - see also xtight.m

  if nargin && isaxis(varargin{1}), l=2;
     ah=varargin{1}; n=numel(ah);
     if n~=1, wbdie('invalid usage (got %d axis handles)',n); end
  else ah=gca; l=1;
  end

  if l>nargin, varargin{1}='-data';
  elseif l<nargin
     if ~helpthis(nargout,varargin{l:end}), wbdie('invalid usage'); end
     return
  end

  switch varargin{l}
    case '-pos',  xx=getXLimData(ah,varargin{l});
    case '-view', xx=getXLimView(ah);
    case '-data', xx=getXLimData(ah);
    otherwise wbdie('invalid usage')
  end

end

% -------------------------------------------------------------------- %
% determine xlim only within ylim view
% skip line handles with viewer than 3 data points
% since quite likely these are just markers

function xx=getXLimView(ah)

  lh=findall(ah,'Type','Line','visible','on');
  m=length(lh); xx=zeros(m,3);

  yl=ylim;

  islog=isequal(get(ah,'XScale'),'log');

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
     xx=xlim(ah);
  end

end

% -------------------------------------------------------------------- %
% NB! on log scale xlim(1) might be =0 even though data won't be plotted !!

function xx=getXLimData(ah,pflag)

  if nargin>1
     if ~isequal(pflag,'-pos'), wbdie('invalid usage'); end
  pflag=1; else pflag=0; end

  lh=findall(ah,'Type','Line','visible','on');
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

