function xl=xtight(varargin)
% function xl=xtight([ah,][fac,][ opts])
%
%    Set x-axis tight to data. Together with ytight(), this 
%    represents a more controlled version to Matlab's 'axis tight'.
%
% Options
%
%   fac      extra zoom (default fac=1, fac>1 zooms out)
%   ah       uses specified axis handle instead of defaul (current axis)
%
%   'view'   zoom to xlim within local xlim view
%   'x1',..  fixed xlim(1)
%   'x2',..  fixed xlim(2)
%
% Wb,2002

  if nargin && isaxis(varargin{1}), l=2;
     ah=varargin{1}; n=numel(ah);
     if n~=1, wbdie('invalid usage (got %d axis handles)',n); end
  else ah=gca; l=1;
  end

  getopt ('init', varargin(l:end));
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
     xl=getxlim(ah,'-view'); dx=diff(xl);
     if dx<=0
        s=sprintf('getxlim() returned [%s]',vec2str(xl));
        if dx<0, wbdie(s); else wblog('WRN',s); end
     else
        xl=xlim__(ah,xl);
     end
  elseif dflag
     xl=xlim__(ah,getxlim(ah,'-data'));
  else
     xopts={'YLim',get(ah,'YLim'),'YLimMode',get(ah,'YLimMode')};
     axis(ah,'tight');
     xl=xlim__(ah,xlim(ah),xopts{:});
  end

  if ~isempty(fac)
     if isequal(get(ah,'XScale'),'linear')
        dx=(fac-1)/2*diff(xl);
        xl=xlim__(ah,[xl(1)-dx, xl(2)+dx]);
     else
        if xl(1)==0, xl=getxlim(ah,'-data'); end
        if all(xl~=0)
           fac=exp((fac-1)/2*diff(log(abs(xl))));
           if fac~=0 && ~any(isinf(xl)) && ~(isnan(fac) || isinf(fac))
           xl=xlim__(ah,[xl(1)/fac, xl(2)*fac]); end
        end
     end
  elseif isequal(get(ah,'YScale'),'log') && xl(1)==0
     xl=xlim__(ah,getxlim(ah,'-pos')); % ,'-data'
  end

  if ~isempty(x1) || ~isempty(x2), xl=xlim(ah);
     if ~isempty(x1), xl(1)=x1; end
     if ~isempty(x2), xl(2)=x2; end
     xl=xlim__(ah,xl);
  end

  if ~isempty(mh), setprops(mh,'-reset'); end
  if nargout, xl=xlim(ah); else clear xl; end

end

function xl=xlim__(ah,xl,varargin)
   if diff(xl)==0, return; end
   xr=round(xl); e=abs(xr-xl)/norm(diff(xl));
   i=find(e<0.025); if ~isempty(i), xl(i)=xr(i); end
   set(ah,'XLim',xl,varargin{:});
end

% -------------------------------------------------------------------- %

