function yl=ytight(varargin)
% function yl=ytight([ah,][fac,][ opts])
%
%    Set y-axis tight to data. Together with xtight(), this
%    represents a more controlled version to Matlab's 'axis tight'.
%
% Options
%
%   fac      extra zoom (default fac=1, fac>1 zooms out)
%   ah       uses specified axis handle instead of defaul (current axis)
%
%   'view'   zoom to ylim within local xlim view
%   'y1',..  fixed ylim(1)
%   'y2',..  fixed ylim(2)
%   'y1>',.. ylim(1) at least greater or equal specified value
%   'y2>',.. ylim(2) at least greater or equal specified value
%   'y1<',.. ylim(1) at least lower or equal specified value
%   'y2<',.. ylim(2) at least lower or equal specified value
%   '-a1'    ytight while keeping DataAspectRatio to 1.
%   '-off'   axis off
%
% Wb,2002

  if nargin && isaxis(varargin{1}), l=2;
     ah=varargin{1}; n=numel(ah);
     if n~=1, wbdie('invalid usage (got %d axis handles)',n); end
  else ah=gca; l=1;
  end

  a1flag=0;

  getopt('init',varargin(l:end));
     y1 = getopt('y1',[]);
     y2 = getopt('y2',[]);
     y1l= getopt('y1<',[]);
     y1g= getopt('y1>',[]);
     y2l= getopt('y2<',[]);
     y2g= getopt('y2>',[]);
     vw = getopt('view');

     if getopt('-a1' ), a1flag=1;
     elseif getopt('-a1L' ), a1flag='L';
     elseif getopt('-a1R' ), a1flag='R'; end
     dx=getopt('dx',0);

     offlag = getopt('-off');
  fac=getopt('get_last',[]);

  mh=[ findall(ah,'tag','xmark'); findall(ah,'tag','ymark') ];
  if ~isempty(mh)
     setprops(mh,'-s','Visible','off');
  end

  if a1flag

     p=get(ah,'Position'); a=p(3)/p(4);
     if ~isempty(fac), fac={fac}; else fac={}; end

     axis(ah,'equal','tight');
        ytight(ah,fac{:});
        xl=xlim(ah); yl=ylim(ah); autoax(ah,'-x');
     set(ah,'PlotBoxAspectRatio',[a 1 1],'DataAspectRatio',[1 1 1],'YLim',yl)
     if a1flag=='L'
        set(ah,'XLim',xl(1)-dx+[0, diff(get(ah,'XLim'))]);
     elseif a1flag=='R'
        set(ah,'XLim',xl(2)+dx-[diff(get(ah,'XLim')), 0]);
     elseif dx
        set(ah,'XLim',dx+[0, diff(get(ah,'XLim'))]);
     end

     if ~isempty(mh), setprops(mh,'-reset'); end
     if ~nargout, clear yl; end
     if offlag, axis(ah,'off'), end
     return
  end

  if ~vw
     xl=xlim(ah); if isequal(get(ah,'XScale'),'log') && xl(1)==0
     xl=getxlim('-data'); end

     axis(ah,'tight');
     if diff(xl)>0
        set(ah,'XLim',xl,'XLimMode', get(ah,'XLimMode'));
     end
     yl=ylim(ah);
  else
     yl=getylim(ah,'-view'); if any(isnan(yl)), return; end
     e=diff(yl);
       if e<0, wbdie('getylim() returned [%s]',vec2str(yl)); end
       if e==0, yl=ylim(ah); return; end
     ylim(ah,yl)
  end

  if ~isempty(fac)

     if isequal(get(ah,'YScale'),'linear')
        dy=(fac-1)/2*diff(yl); set(ah,'YLim', [yl(1)-dy, yl(2)+dy ]);
     else
        if yl(1)==0, yl=getylim(ah,'-data'); end
        if all(yl~=0)
           fac=exp((fac-1)/2*diff(log(abs(yl))));
           if fac~=0 && ~any(isinf(yl)) && ~(isnan(fac) || isinf(fac))
           set(ah,'YLim', [yl(1)/fac, yl(2)*fac]); end
        end
     end

  elseif isequal(get(ah,'YScale'),'log') && yl(1)==0
     set(ah,'YLim',getylim(ah,'-data'));
  end

  yl=ylim(ah);
     if ~isempty(y1), yl(1)=y1;
     else
        if ~isempty(y1g), if y1g>yl(1), yl(1)=y1g; end, end
        if ~isempty(y1l), if y1l<yl(1), yl(1)=y1l; end, end
     end

     if ~isempty(y2), yl(2)=y2;
     else
        if ~isempty(y2g), if y2g>yl(2), yl(2)=y2g; end, end
        if ~isempty(y2l), if y2l<yl(2), yl(2)=y2l; end, end
     end
  ylim(ah,yl);

  if ~isempty(mh), setprops(mh,'-reset'); end
  if nargout, yl=ylim(ah); else clear yl; end

end

% -------------------------------------------------------------------- %

