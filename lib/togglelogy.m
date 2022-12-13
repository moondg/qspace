function togglelogy(varargin)
% function togglelogy([ah,opts])
%
%    toggle y-scale between log|y| and lin(y)
%    data for negative y values is blurred (shown in lighter color)
%
% Options
%
%    ~b    do not blur data for negative y values
%
% Wb,Oct28,05 ; Wb,Dec05,19

  getopt('init',varargin);
     noblur=getopt('~b');
  ah=getopt('get_last',[]); if isempty(ah), ah=gca; end

  tag=mfilename;
  ntag=[tag ':bndata'];

  ish=ishold;
  makelog=~isequal(get(ah,'YScale'),'log');

  if makelog
      au=get(ah,'UserData'); if ~isstruct(au), au=struct; end
      au.xlm=get(ah,'XLim');
      au.ylm=get(ah,'YLim'); set(ah,'UserData',au);

      hh=findall_lines(ah);
      for h=hh
          u=get(h,'UserData'); if ~isstruct(u), u=struct; end

          y=get(h,'YData');
          in=find(y<=0); if isempty(in)
             setfield(u,tag,0); set(h,'UserData',u);
             continue
          end

          u.yd=y;
          if noblur, y=abs(y);
          else
             xd=get(h,'XData'); o=lineopts(h); hold on
             y(find(y>=0))=nan;
                u.h2=plot(xd,-y,o{:},'tag',ntag);
                mv2back(blurl(u.h2,'cfac',0.66,'LineW',1));
             y=u.yd; y(in)=nan;
          end
          set(h,'YData',y,'UserData',u);
      end

      set(ah,'YLimMode','auto','YScale','log'); axis tight
      yl=ylim; yl(2)=max(abs(au.ylm)); ylim(yl); xlim(au.xlm);

      if ~ish, hold off; end
  else
      set(ah,'YScale','lin');
      au=get(ah,'UserData'); e=0;

      if all(isfield(au,{'xlm','ylm'}))
         set(ah,'XLim', au.xlm,'YLim', au.ylm);
         delete(findall(ah,'tag',ntag));

         hh=findall_lines(ah); nh=numel(hh);
         for h=hh
             u=get(h,'UserData');
             if isfield(u,'yd'), set(h,'YData',u.yd);
             elseif ~isfield(u,tag) || getfield(u,tag), e=e+1; end
         end
         if e, wblog('ERR missing lin-scale UserData (e=%g/%g)',e,nh); end
     end
  end

  if exist('ymark','file'), ymark -fix; end

end

