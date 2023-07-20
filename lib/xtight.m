function xl=xtight(varargin)
% xtight - set x axis tight to data
% usage: xtight([fac, opts])
%
%    fac = extra zoom (default fac=1, fac>1 introduces white margin)
%
% Options
%
%   'view'   zoom to xlim within local ylim view
%   'x1',..  fixed xlim(1)
%   'x2',..  fixed xlim(2)
%
% Wb,2002

  if nargin && isaxis(varargin{1})
       ah=varargin{1}; varargin(1)=[]; setax(ah);
  else ah=gca; end

  getopt ('init', varargin);
     x1 = getopt('x1',[]);
     x2 = getopt('x2',[]);
     vw = getopt('view' );
     dflag= getopt('-d' );
  fac=getopt('get_last',[]);

  if ~isnumeric(fac), wbdie('invalid usage (fac)'); end

  mh=[ findall(ah,'tag','xmark'); findall(ah,'tag','ymark') ];
  if ~isempty(mh)
     setprops(mh,'-s','Visible','off');
  end

  if vw
     xl=getxlim('-view'); dx=diff(xl);
     if dx<=0
        s=sprintf('getxlim() returned [%s]',vec2str(xl));
        if dx<0, wbdie(s); else wblog('WRN',s); end
     else
        xl=xlim__(ah,xl);
     end
  elseif dflag
     xl=xlim__(ah,getxlim('-data'));
  else
     xopts={'YLim',get(ah,'YLim'),'YLimMode',get(ah,'YLimMode')};
     axis tight;
     xl=xlim__(ah,xlim,xopts{:});
  end

  if ~isempty(fac)
     if isequal(get(ah,'XScale'),'linear')
        dx=(fac-1)/2*diff(xl);
        xl=xlim__(ah,[xl(1)-dx, xl(2)+dx]);
     else
        if xl(1)==0, xl=getxlim('-data'); end
        if all(xl~=0)
           fac=exp((fac-1)/2*diff(log(abs(xl))));
           if fac~=0 && ~any(isinf(xl)) && ~(isnan(fac) || isinf(fac))
           xl=xlim__(ah,[xl(1)/fac, xl(2)*fac]); end
        end
     end
  elseif isequal(get(ah,'YScale'),'log') && xl(1)==0
     xl=xlim__(ah,getxlim('-data','-pos'));
  end

  if ~isempty(x1) || ~isempty(x2), xl=xlim;
     if ~isempty(x1), xl(1)=x1; end
     if ~isempty(x2), xl(2)=x2; end
     xl=xlim__(ah,xl);
  end

  if ~isempty(mh), setprops(mh,'-reset'); end
  if nargout, xl=xlim; else clear xl; end

end

function xl=xlim__(ah,xl,varargin)
   xr=round(xl); e=abs(xr-xl)/norm(diff(xl));
   i=find(e<0.025); if ~isempty(i), xl(i)=xr(i); end
   set(ah,'XLim',xl,varargin{:});
end

% -------------------------------------------------------------------- %

