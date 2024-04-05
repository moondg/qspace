function xx=getylim(varargin)
% Function xx=getylim([ah,][opts])
% Options
%
%    '-view'   determine YLim fully shows all y-data in given view
%    '-data'   default: ylim is full y-data range
%
% Wb,2002; Wb,Jan10,08 - see also ytight.m

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
    case '-view', xx=getYLimView(ah);
    case '-data', xx=getYLimData(ah);
    otherwise wbdie('invalid usage')
  end

end

% -------------------------------------------------------------------- %
% determine ylim only within xlim view
% skip line handles with viewer than 3 data points
% since quite likely these are just markers

function yy=getYLimView(ah)

  lh=findall(ah,'Type','Line','visible','on');
  m=length(lh); if ~m, yy=[]; return; end

  yy=nan(m,3); xl=xlim(ah);

  for i=1:m
      xd=get(lh(i),'XData'); ii=find(xd>=xl(1) & xd<=xl(2));
      if isempty(ii), continue; end

      yd=get(lh(i),'YData'); yd=yd(ii);
      yy(i,:)=[min(yd), max(yd), length(xd)];
  end

  i=find(yy(:,3)>2); if isempty(i), i=1:size(yy,1); end
  yy=[ min(yy(i,1)), max(yy(i,2)) ];

  if yy(1)==1E99
     wblog('WRN','failed to determine YLim for data within view');
     yy=ylim(ah);
  end
end

% -------------------------------------------------------------------- %
% NB! on log scale ylim(1) might be =0 even though data won't be plotted !!

function yy=getYLimData(ah)

  lh=findall(ah,'Type','Line','visible','on');
  m=length(lh); yy=zeros(m,2);

  for i=1:m
      yd=get(lh(i),'YData'); yd=yd(find(yd>0));
      if ~isempty(yd)
           yy(i,:)=[min(yd), max(yd)];
      else yy(i,:)=[+inf -inf];
      end
  end

  yy=[ min(yy(:,1)), max(yy(:,2)) ];

end

% -------------------------------------------------------------------- %

