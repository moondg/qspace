function cl=clim(varargin)
% Function cl=clim([ah,] OPTS)
%
%   similar behavior as [xyz]lim but for clim
%   including: clim auto
%
% Additional options
%
%   'cfac',r  scale CLim to data range; r=1 corresponds
%             to exact data range, while r<1 truncates.
%   'c[12]',. fixed value for CLim([12])
%   
% Wb,Aug13,08

  if nargin && isaxis(varargin{1})
     ah=varargin{1}; varargin(1)=[];
  else ah=gca; end

  getopt('init',varargin);
     aflag=getopt('auto');
     cfac =getopt('cfac',[]);
     c1   =getopt('c1',[]);
     c2   =getopt('c2',[]);
  varargin=getopt('get_remaining'); narg=length(varargin);

  if narg>1 || narg && ~isnumeric(varargin{1})
     eval(['help ' mfilename]);
     if nargin || nargout, wbdie('invalid usage'), end, return
  end

  if narg, cl=varargin{1}; else cl=get(ah,'CLim'); end
  if aflag || nargin==0
     if nargin && ~nargout, clear cl; end
     return
  end

  if aflag, set(ah,'CLimMode','auto'); end
  if ~isempty(cfac)
     sh=findall(ah,'type','surface'); nh=length(sh); I=cell(nh,1);
     for i=1:nh, I{i}=reshape(get(sh,'ZData'),[],1); end

     d=cat(1,I{:}); if ~isempty(d), d=sort(d); n=numel(d);
       if isscalar(cfac), cfac=[max([0,0.5*(1-cfac)]), 0.5*(1+cfac)]; end;
       cl=d(ceil(n*cfac));
     end
  end

  if ~isempty(c1), cl(1)=c1; end
  if ~isempty(c2), cl(2)=c2; end

  set(ah,'CLim',cl); if ~nargout, clear cl; end

end

