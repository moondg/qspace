function [h,th]=ymark(varargin)
% usage: ymark(ypos [,opts])
% mark y position in current figure
%
% Options:
%
%   'istr',str  info string positioned by default at relative xpos=0.12
%   'istr',{str, [xpos,][text options]}
%
%   '-rm'   remove all existing ymarks
%   '-fix'  redraw ymarks so they stretch over full xlim-range
%   '-bg'   move ymark into background
%   'xl',.. limits of horizontal marker in normalized units
%   'zd',.. specify ZData
%   '-0'    show axis at y=0
%
%   remainder of options are used as options to ymark line to be drawn.
%
% See also xmark, allmarks, ypatch.

  getopt('init',varargin);
     istr =getopt('istr','');
     rmall=getopt('-rm'  );
     fxall=getopt('-fix');
     bflag=getopt('-bg');
     fgflg=getopt('-fg');
     xl1  =getopt('xl',[]);
     zd   =getopt('zd',[]);
     y0   =getopt('-0');
  varargin=getopt('get_remaining'); narg=numel(varargin);

  if rmall
     delete(findall(gca,'tag','ymark'));
     if ~narg, return; end
  end

  xl=xlim; gotlog=isequal(get(gca,'XScale'),'log');
  if gotlog && xl(1)==0
     t=get(gca,'XTick'); if length(t>1), xl(1)=t(1);
     else wblog('WRN','xlim returns 0 as left boundary'); end
  end

  if ~isempty(xl1)
     if numel(xl1)~=2, wbdie('invalid xl(im) values'); end
     if gotlog, xl=log(xl); end
     xl=xl(1)+xl1*diff(xl);
     if gotlog, xl=exp(xl); end
  end

  if fxall
     for m=findall(gca,'type','line','tag','ymark')', set(m,'XData',xl); end
     if ~narg, return; end
  end

  if y0, y=ylim; if diff(sign(y))==2
     if isempty(varargin)
          varargin={0};
     elseif isnumeric(varargin{1})
          varargin{1}=[varargin{1}, 0];
     else varargin=[0, varargin];
     end
     if numel(varargin{1})==1, bflag=1; end
     if numel(varargin)<2
        varargin=[varargin,{'Color',[.9 .9 .9], 'LineWidth',0.75}];
     end
  end, end

  if isempty(varargin), return; end
  if isnumeric(varargin{1})
     y=varargin{1}; varargin(1)=[];
  end

  hold on

  if mod(numel(varargin),2)
       fmt=varargin(1); varargin(1)=[];
  else fmt={'Color', [1  1 .7],'LineW',4}; bflag=1; end

  if fgflg, bflag=0; end

  h=plot(xl, [y(:) y(:)]', fmt{:},'tag','ymark');
  if ~isempty(zd)
     set(h,'ZData',zd([1 1]));
  elseif bflag
     mv2back(h);
     set(h,'ZData',repmat(-0.1,size(get(h,'XData'))));
  end

  if ~isempty(varargin)
  set(h,varargin{:}); end

  if ~isempty(istr)
     xt=0.12; n=numel(y); th=zeros(n,1);
     fw={}; % {'FontWeight','bold'};

     if iscell(istr), topt=istr(2:end); istr=istr{1};
        if ~isempty(topt) && isnumeric(topt{1})
           xt=topt{1}; topt(1)=[];
        end
     else topt={}; end

     c=get(h(1),'Color');
     if norm(c)>0.8, c=0.8*c/norm(c); end

     for i=1:n
        th(i)=text(1,y(i),istr,'Color',c,fw{:},'tag','ymark',topt{:});
        set(th(i),'Units','norm'); p=get(th(i),'Pos'); p(1)=xt;
        set(th(i),'Pos',p,'Units','data');
        if ~isempty(zd)
           p=get(th(i),'Pos'); p(3)=zd(1); set(th(i),'Pos',p);
        end
     end
  end

  if ~nargout, clear h; end

end

