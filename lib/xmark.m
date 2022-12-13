function [h,th]=xmark(varargin)
% function [h,th]=xmark(xpos [,opts]);
%
%    Mark specified x-positions in current axis.
%
% Options
%
%   'istr',..   display istr with last xmark line
%   'itop'      display istr on upper end of xmark line (default: bottom)
%   'iTOP'      display istr on top (outside axis set)
%   'yl',[]     ylim of xmark lines in normalized units
%   '-fg'       keep xmark lines front (default: moved to background)
%   '-bg'       move ymark into background
%   '-rm'       remove all marks
%   'zd',..     specify ZData
%   '-0'        show axis at x=0
%
% Example: put 'istr' at ypos=0.8 (normalized) in vertical orientation
%
%   'istr',{'t_{Kondo}',0.8,'Rotation',90})
%
% See also allmarks, xpatch.
% Wb,Mar07

  bflag=1;

  getopt('INIT',varargin);
     istr =getopt('istr','');
     itop =getopt('itop',''); if ~isempty(itop), istr=itop; itop=1; end
     iTOP =getopt('iTOP',''); if ~isempty(iTOP), istr=iTOP; itop=2; end
     if  getopt({'-fg','top'}), bflag=0;
     else  getopt('-bg'); end 
     yl   =getopt('yl',[]);
     zd   =getopt('zd',[]);
     rmall=getopt('-rm');
     gflag=getopt('-g' );
     fxall=getopt('-fix');
     r90  =getopt('-90' );
     x0   =getopt('-0');
  varargin=getopt('get_remaining');

  if rmall
     delete(findall(gca,'tag','xmark')); % 'type','line'
     if nargin<2, return; end
  end

  dm=0.005;

  if fxall
     yl=ylim; s=get(gca,'YScale');
     if isequal(s,'log'), if yl(1)==0, yl(1)=1E-16; end
          dy=(yl(2)/yl(1))^dm; yl=[yl(1)*dy, yl(2)/dy];
     else dy=dm*diff(yl); yl=[yl(1)+dy, yl(2)-dy]; end

     hh=findall(gca,'type','line','tag','xmark')';
     for m=hh, set(m,'YData',yl); end
     if nargin<2, return; end
  end

  if x0, x=xlim; if diff(sign(x))==2
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
     x=varargin{1};
     varargin(1)=[];
  else x=[]; end

  hold on

  if mod(length(varargin),2)
       fmt=varargin(1); varargin(1)=[];
  else fmt={'Color', [1  1 .7],'LineW',4}; end

  if isempty(yl)
     yl=ylim; s=get(gca,'YScale');
     if isequal(s,'log')
          dy=(yl(2)/yl(1))^dm; yl=[yl(1)*dy, yl(2)/dy];
     else dy=dm*diff(yl); yl=[yl(1)+dy, yl(2)-dy]; end
  else l=ylim;
     if isequal(get(gca,'YScale'),'log') && all(l>0)
          yl=l(1)*(l(2)/l(1)).^yl;
     else yl=l(1)+diff(l).*yl; end
  end

  if isempty(x), return; end

  if any(isnan(yl))
     yl=getylim('-data');
  end

  h=plot([x(:) x(:)]', yl,fmt{:},'tag','xmark');
  if ~isempty(zd)
     set(h,'ZData',zd([1 1]));
  elseif bflag, mv2back(h); end

  if gflag
  set(h,'Color',[1 1 1]*.9, 'LineW',0.5); end

  if ~isempty(varargin)
  set(h,varargin{:}); end

  if ~isempty(istr)
     c=0.7*get(h(1),'Color');
     fw={}; % {'FontWeight','bold'};
     bw={}; yt=0.12;
     if r90, bw={'Rotation',90,'BackgroundC','w','Margin',1E-3}; end

     if iscell(istr), topt=istr(2:end); istr=istr{1};
        if length(topt) && isnumeric(topt{1})
        yt=topt{1}; topt=topt(2:end); end
     else topt={}; end

     switch itop
        case 1, yt=0.90; 
        case 2, yt=1.15; c='k'; fw={};
     end

     n=numel(x); th=[];

     xl=xlim; if x(end)>=xl(1) && x(end)<=xl(2)
        th=text(x(end),1,istr,'Color',c,fw{:},bw{:},...
          'tag','xmark','UserD',yt,'HorizontalAl','center',topt{:});
        set(th,'Units','norm'); p=get(th,'Pos'); p(2)=yt;
        set(th,'Pos',p); set(th,'Units','data');
        if ~isempty(zd)
           p=get(th,'Pos'); p(3)=zd(1); set(th,'Pos',p);
        end
     end
  end

  if ~nargout, clear h; end

end

