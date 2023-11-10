function th=show_deg(Ek,varargin)
% Function th=show_deg(Ek [,ah,opts])
%
%    show degeneracies of line plot (e.g. NRG flow diagram).
%
% Options
%
%    'dz',..    internal multiplet dimension (one row for each value Ek)
%    'qq',..    corresponding symmetry labels (one row for each value Ek)
%
%    'fs',..    fontsize of labels (8)
%    'tag',...  tag to apply to text labels (':deg:')
%    'fs',..    fontsize of labels (8)
%    'eps',...  consider levels degenerate within eps (1E-3)
%    'x0',...   x-position of labels (if 0<x0<1, assume relative x-axis position)
%    'xt',...   same as x0, but always in absolute units
%    'dx',...   x-offset if non-degenerate y-values by less than dy ([dxmin,dxmax])
%    'dy',...   threshold to use dx (see dx)
%    'n',...    maximum number of values to consider in sorted Ek
%    'y2',...   ignore values Ek>y2 (default: 0.6 of ylim range)
%   
%    'nmax', .. do not show if number of closeby levels exceeds nmax (Inf)
%    'lmax', .. show at most the lowest lmax levels for each symmetry (Inf)
%   
%     --no1s    skip labels that show `1'
%
% see also nrg_plot.m
% Wb,Mar18,08

  if nargin>1 && isaxis(varargin{1})
     setax(varargin{1}); varargin=varargin(2:end);
  end

  fs=8;
  getopt('init',varargin);

     dz   = getopt('dz',[]);
     qq   = getopt('qq',[]);
     vflag= getopt('-v');

     tag  = getopt('tag',':deg:');
     fs   = getopt('fs',fs);
     xt   = getopt('xt',[]); if isempty(xt)
     x0   = getopt('x0',0.8); end
     dx   = getopt('dx',[]);
     dy   = getopt('dy',0.25);
     nmax = getopt('nmax',Inf);
     lmax = getopt('lmax',Inf);
     x2   = getopt('x2',[]);
     y2   = getopt('y2',[]);
     eps  = getopt('eps',1E-3);
     n    = getopt('n',128);
     no1s = getopt('--no1s');
  varargin= getopt('get_remaining'); narg=length(varargin);

  if numel(y2)==1 && y2>0, dy=dy*y2; end

  if isempty(qq)
       va={'VerticalAl','bottom','Margin',1E-3}; % ,'BackgroundC','w'
  else va={'VerticalAl','bottom'}; end
  topt={'HorizontalAl','center',va{:},'tag',tag,'FontSize',fs,varargin{:}};

  yl=ylim; if isempty(y2), y2=yl(1)+0.6*diff(yl); end
  n=numel(Ek);

  if ~isvector(Ek)
     wbdie('invalid usage (expecting vector for Ek)'); end
  if ~isempty(dz) && size(dz,1)~=n, whos dz Ek
     wbdie('invalid usage (size mismatch of dz)'); end
  if ~isempty(qq) && size(qq,1)~=n, whos qq Ek
     wbdie('invalid usage (size mismatch of qq)'); end

  [Ek,is]=sort(reshape(Ek,[],1));
  Ek(find(isnan(Ek)))=Inf;

  if ~isempty(dz), dz=dz(is,:); end
  if ~isempty(qq), qq=qq(is,:); end

  Ek(end+1)=max(y2,Ek(end))+10*eps;

  I=[1; 1+find(diff(Ek)>=eps) ];
  i=find(Ek(I)>y2); I(i(2:end))=[];

  if Ek(I(end))>y2+20*eps, I(end)=[]; end

  xl=xlim;
  if ~isempty(xt), x0=xt;
  elseif x0<=1, x0=xl(1)+x0*diff(xl); end

  if isempty(dx), dx=diff(xl)/50; end
  m=length(I); xi=repmat(x0,m,1);

  Ex=Ek(I); i1=1; n2=numel(Ex); nx=numel(dx); count=0;
  while i1<m
     i2=i1+find(Ex(i1+1:n2)>(Ex(i1)+dy),1);
     if isempty(i2)
        if i1<n2, i2=n2; else break; end
     end

     count=count+1;
     if count>lmax, xi(i1:end)=Inf; break; end

     if i2>i1+1, q=(i2-i1-1)/2; q=(-q:q);
        if nx==1, q=q*dx;
        else
           nq=numel(q); if nq>nmax, q(:)=Inf;
           else
              if nq==2, q=q*dx(2);
              elseif nq==3, q=q*mean(dx(1:2));
              else q=q*dx(1); end
           end
        end
        xi(i1:(i2-1))=x0+q;
     end
     i1=i2;
  end

  de=0.1*mean(diff(Ex));

  for k=2:m
     if ~isempty(x2) && xi(k-1)>x2, continue; end

     i=I(k-1):(I(k)-1);
     d=numel(i); if d==1 && no1s, continue; end
     if isempty(dz)
        dstr=sprintf('%g',d);
     else
        g=prod(dz(i,:),2);
        if numel(g)>1 && vflag
           s=sprintf('%+g',g);
           dstr=sprintf('%s=%g',s(2:end),sum(g));
        else dstr=sprintf('%g',sum(g)); end
     end
     th(k-1,1)=text(xi(k-1),mean(Ek(i)),1, dstr,topt{:},'VerticalA','middle');

     if ~isempty(qq)
        if isempty(dz)
           th(k-1,2)=text(xi(k-1),mean(Ek(i)), ...
              mat2str2(qq(i,:),'fmt','%2g'), ...
              topt{:},'VerticalA','top','FontSize',fs-2 ...
           );
        else
           fmt=[ '[ ' repmat('%2g ',1,size(qq,2)) '] (%g)\n' ]; g=prod(dz(i,:),2);
           th(k-1,2)=text(xi(k-1),mean(Ek(i)), ...
              sprintf(fmt,[qq(i,:),g]'),...
              topt{:},'VerticalA','top','FontSize',fs-2 ...
           );
        end
     end
  end

  if ~nargout, clear th; end

end

