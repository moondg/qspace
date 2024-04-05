function lh=scalex(varargin)
% function lh=scalex([ah,] sc [,opts])
%
%    change current units of x-axis x -> x*sc 
%    for axis ah (default: gca)
%
% Options:
%
%    '-p'        interpret sc as polynomial coefficients
%    'xstr',..   new xlabel to use
%
%    'glines'    draw vertical grid lines at positions indicated (grl in units of sc)
%    'ulabel',.. unit label
%    'ilabel',.. inverse unit label (here: divide by sc!)
%
% Wb,Oct24,05

  if nargin && isaxis(varargin{1}), l=2;
     ah=varargin{1}; n=numel(ah);
     if n~=1, wbdie('invalid usage (got %d axis handles)',n); end
  else ah=gca; l=1;
  end

  if nargin>=l, sc=varargin{2}; l=l+1;
  else wbdie('invalid usage'); end

  pflag=numel(sc)>1;
  ulb=''; grl=[]; iflag=0;

  getopt('init',varargin(l:end))
     if pflag
        ulb = getopt('xstr',[]);
        if numel(sc)<2, wbdie('got invalid (constant) polynomial!'); end
     else
        if numel(sc)~=1, wbdie('got invalid scale factor!'); end
        grl=getopt('glines',grl);
        ulb=getopt('ulabel',ulb);
        if isempty(ulb)
           ulb=getopt('ilabel',ulb);
           if ~isempty(ulb), iflag=1; end
        end
     end
  varargin=getopt('get_remaining');

  if ~isempty(varargin)
     if isempty(ulb) && length(varargin)==1 && ischar(varargin{1})
          ulb=varargin{1};
     else wbdie('invalid usage'); end
  end

  if iflag, sc=1/sc; end

  if ~isempty(ulb)
     th=get(ah,'XLabel'); s=get(th,'string');
     if pflag, s=ulb; else
        if iflag
           s=[s '/' ulb];
        else
           s=[s '{\cdot}' ulb];
        end
     end
     set(th,'String',s)
  end

  set(ah,'XLim',scale_data(get(ah,'XLim'),sc));

  for lh=[ findall(ah,'type','line'); findall(ah,'type','patch') ]'
     set(lh,'XData', scale_data(get(lh,'XData'),sc));
  end

  for h=[ findall(ah,'type','text') ]'
     if ~isequal(get(h,'units'),'data') || islabel(h), continue, end
     p=get(h,'Position'); p(1)=scale_data(p(1),sc);
     set(h,'Position',p);
  end

  if ~isempty(grl)
     hold on
     n=length(grl); lh=zeros(length(n));
     for i=1:n
         lh(i)=plot([grl(i), grl(i)], ylim, 'k--');
     end
     xl=[min(grl), max(grl)]; xl = xl+0.2*[-1 1].*diff(xl); xlim(xl);
  end

  set(ah,'XTickMode','auto');

  if ~nargout, clear lh; end

end

% -------------------------------------------------------------------- %
% Wb,Mar14,15

function x=scale_data(x,sc)
   if numel(sc)==1
        x=sc*x;
   else x=polyval(sc,x);
   end
end

% -------------------------------------------------------------------- %

