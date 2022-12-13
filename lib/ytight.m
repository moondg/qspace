function yl=ytight(varargin)
% ytight - set y axis tight to data
% usage: ytight([fac, opts])
%
%    fac = extra zoom (default fac=1, fac>1 introduces white margin)
%
% Options
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

  a1flag=0;

  getopt('init',varargin);
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

  ah=gca;

  mh=[ findall(ah,'tag','xmark'); findall(ah,'tag','ymark') ];
  if ~isempty(mh)
     setprops(mh,'-s','Visible','off');
  end

  if a1flag

     p=get(gca,'Position'); a=p(3)/p(4);
     if ~isempty(fac), fac={fac}; else fac={}; end

     axis equal tight; ytight(fac{:}); xl=xlim; yl=ylim; autoax -x
     set(ah,'PlotBoxAspectRatio',[a 1 1],'DataAspectRatio',[1 1 1],'YLim',yl)
     if a1flag=='L'
        set(gca,'XLim',xl(1)-dx+[0, diff(get(gca,'XLim'))]);
     elseif a1flag=='R'
        set(gca,'XLim',xl(2)+dx-[diff(get(gca,'XLim')), 0]);
     elseif dx
        set(gca,'XLim',dx+[0, diff(get(gca,'XLim'))]);
     end

     if ~isempty(mh), setprops(mh,'-reset'); end
     if ~nargout, clear yl; end
     if offlag, axis off, end
     return
  end

  if ~vw
     xl=xlim; if isequal(get(ah,'XScale'),'log') && xl(1)==0
     xl=getxlim('-data'); end

     xopts={'XLim',xl,'XLimMode', get(ah,'XLimMode')};
     axis tight; set(ah,xopts{:}); yl=ylim;
  else
     yl=getylim('-view'); if any(isnan(yl)), return; end
     e=diff(yl);
       if e<0, wberr('getylim() returned [%s]',vec2str(yl)); end
       if e==0, yl=ylim; return; end
     ylim(yl)
  end

  if ~isempty(fac)

     if isequal(get(ah,'YScale'),'linear')
        dy=(fac-1)/2*diff(yl); set(ah,'YLim', [yl(1)-dy, yl(2)+dy ]);
     else
        if yl(1)==0, yl=getylim('-data'); end
        if all(yl~=0)
           fac=exp((fac-1)/2*diff(log(abs(yl))));
           if fac~=0 && ~any(isinf(yl)) && ~(isnan(fac) || isinf(fac))
           set(ah,'YLim', [yl(1)/fac, yl(2)*fac]); end
        end
     end

  elseif isequal(get(ah,'YScale'),'log') && yl(1)==0
     set(ah,'YLim',getylim('-data'));
  end

  yl=ylim;
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
  ylim(yl);

  if ~isempty(mh), setprops(mh,'-reset'); end
  if nargout, yl=ylim; else clear yl; end

end

% -------------------------------------------------------------------- %

