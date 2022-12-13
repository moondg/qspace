function togglelogx(varargin)
% function togglelogx([ah,opts])
%
%    toggle x-scale between log|x| and lin(x)
%    data for negative x values is blurred (shown in lighter color)
%
% Options
%
%    ~b    do not blur data for negative x values
%
% Wb,Oct28,05 ; Wb,Dec05,19

  getopt('init',varargin);
     noblur=getopt('~b');
  ah=getopt('get_last',[]); if isempty(ah), ah=gca; end

  tag=mfilename;
  ntag=[tag ':bndata'];

  ish=ishold;
  makelog=~isequal(get(ah,'XScale'),'log');

  if makelog
      au=get(ah,'UserData'); if ~isstruct(au), au=struct; end
      au.xlm=get(ah,'XLim');
      au.ylm=get(ah,'YLim'); set(ah,'UserData',au);

      for h=findall_lines(ah)
          u=get(h,'UserData'); if ~isstruct(u), u=struct; end

          x=get(h,'XData');
          in=find(x<=0); if isempty(in)
             setfield(u,tag,0); set(h,'UserData',u);
             continue
          end

          u.xd=x;
          if noblur, x=abs(x);
          else
             yd=get(h,'YData'); o=lineopts(h); hold on
             x(find(x>=0))=nan;
                u.h2=plot(-x,yd,o{:},'tag',ntag);
                mv2back(blurl(u.h2,'cfac',0.66,'LineW',1));
             x=u.xd; x(in)=nan;
          end
          set(h,'XData',x,'UserData',u);
      end

      set(ah,'XLimMode','auto','XScale','log'); axis tight
      xl=xlim; xl(2)=max(abs(au.xlm)); xlim(xl); ylim(au.ylm);

      if ~ish, hold off; end
  else
      set(ah,'XScale','lin');
      au=get(ah,'UserData'); e=0;

      if all(isfield(au,{'xlm','ylm'}))
         set(ah,'XLim', au.xlm, 'YLim', au.ylm);
         delete(findall(ah,'tag',ntag));

         hh=findall_lines(ah); nh=numel(hh);
         for h=hh
             u=get(h,'UserData');
             if isfield(u,'xd'), set(h,'XData',u.xd);
             elseif ~isfield(u,tag) || getfield(u,tag), e=e+1; end
         end
         if e, wblog('ERR missing lin-scale UserData (e=%g/%g)',e,nh); end
     end
  end

  if exist('ymark','file'), ymark -fix; end

end

