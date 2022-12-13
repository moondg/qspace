function plotlevels(X,Y,varargin)
% function plotlevels(x,y [,opts])
%
% Options
%
%    'dx'        width of level indicator
%    'yl'        ylim ([-Inf Inf])
%    'y2'        show degeneracies only up to y<y2 (default: yl(2))
%    '-u'        show degeneracies only for degenerate spaces
%    '-S'        show total degeneracy at given energy to the right of the plot
%    'topt',{..} options for text on degeneracies
%
% Remaining options are used as line options.
% Wb,Jun01,11

% tags: plotspectrum, plotspectra, plotdeg

  if nargin<2
     helpthis, if nargin || nargout
     wberr('invalid usage'), end, return
  elseif ~isvector(X) || ~isvector(Y)
     wberr('invalid usage'),
  end

  getopt('init',varargin);
     yl=getopt('yl',[]); if ~isempty(yl), y2=yl(2); else y2=Inf; end
     y2=getopt('y2',y2);
     Sflag=getopt('-S');
     uflag=getopt('-u');
     dx=getopt('dx',0.3);
     topt=getopt('topt',{});
  varargin=getopt('get_remaining'); narg=length(varargin);

  if narg && isempty(topt) && iscell(varargin{end})
     topt=varargin{end}; varargin(end)=[];
  end

  if isempty(varargin)
       lopt={'Color',[1 1 1]*0.75};
  else lopt=varargin; end

  if isempty(topt), topt={'FontSize',8}; end
  topt={'HorizontalAl','center','VerticalAl','middle',topt{:}};

  [d1,I1,D1]=uniquerows(X);

  X=X(:); Y=Y(:); dx=0.5*abs(dx);

  for i=1:numel(I1)
     [d2,I2,D2]=uniquerows(Y(I1{i}),'-1');

        j=I1{i}(I2); z=nan(numel(j),1);
        x=X(j); x=reshape([x-dx,x+dx,z]',[],1);
        y=Y(j); y=reshape([y,y,z]',[],1);

     plot(x,y,lopt{:}); hold on

     for j=1:numel(D2)
        if d2(j)<y2 && (~uflag || D2(j)>1)
           text(d1(i),d2(j),sprintf('%g',D2(j)),topt{:});
        end
     end
  end

  if ~isempty(yl), ylim(yl); end

  if Sflag
     x=xlim; xlim(x); x=x(2)+0.02*diff(x); o={'HorizontalAl','left'};
     [y,I,D]=uniquerows(Y);
     for i=1:numel(I)
        if y(i)<y2
           text(x,y(i),sprintf('%g',D(i)),o{:});
        end
     end
  end

end

