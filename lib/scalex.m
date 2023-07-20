function lh=scalex(sc,varargin)
% function lh=scalex(sc [,opts])
%
%    change current units for x-axis: x -> x*sc
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

  pflag=numel(sc)>1;
  ulb=''; grl=[]; iflag=0;

  getopt ('init', varargin);
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
     th=get(gca,'XLabel'); s=get(th,'string');
     if pflag, s=ulb; else
        if iflag
           s=[s '/' ulb];
        else
           s=[s '{\cdot}' ulb];
        end
     end
     set(th,'String',s)
  end

  set(gca,'XLim', scale_data(get(gca,'XLim'),sc));

  for lh=[ findall(gca,'type','line'); findall(gca,'type','patch') ]'
     set(lh,'XData', scale_data(get(lh,'XData'),sc));
  end

  for h=[ findall(gca,'type','text') ]'
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

  set(gca,'XTickMode','auto');

  if ~nargout, clear lh; end

end

% -------------------------------------------------------------------- %
function x=scale_data(x,sc)
   if numel(sc)==1
        x=sc*x;
   else x=polyval(sc,x);
   end
end

% -------------------------------------------------------------------- %

