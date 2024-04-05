function h=inset(varargin)
% Function h=inset([ah,] pos [,opts])
%
%    Position inset to current axis (or to ah's given)
%
% Options
%
%    'scale',fac
%    'dx',..  (in normalized units in current axis set)
%    'ah',..  duplicate data in specified axis handle in inset
%
%     remaining options are applied to axis object.
%
% Wb,Aug15,08

  if nargin && isaxis(varargin{1})
  h0=varargin{1}; varargin(1)=[];
  else h0=gca; end

  getopt('init',varargin);
     rmflag=getopt('-rm');
     dx    =getopt('dx',[]);
     ah    =getopt('ah',[]);
     scale =getopt('scale',1);
  varargin=getopt('get_remaining'); narg=length(varargin);

  tag='inset';
  if rmflag
     delete(findall(h0,'type','text','tag',tag));
     if ~narg, return; end
  end

  if narg<1
     eval(['help ' mfilename]);
     if nargin || nargout, wbdie('invalid usage'), end, return
  end

  pos=postrans(varargin{1});
  if ~isempty(dx)
     if numel(dx)~=2, wbdie('invalid dx'); end
     pos=pos+reshape(dx,1,[]);
  end

  h=h0;
  for i=1:numel(h0)
     h(i)=inset_1(h0(i),pos,scale);

     if ~isempty(ah), j=min(numel(ah),i);
        h(i)=repax(h(i),ah(j));
     end
  end

  set(h,'tag',tag,varargin{2:end});

end

% -------------------------------------------------------------------- %

function h=inset_1(ah,pos,scale)
   if numel(scale)==1, scale=repmat(scale,1,2); end

   p0=get(ah,'Position');
   p(3:4)=scale.*p0(3:4)./[4, 2.75];

   x=pos(1); y=pos(2); if nargin<4, m=0.1; end

   if x<0.25, f=0; elseif x>0.75, f=1; else f=0.5; end
   p(1)=p0(1)+x*p0(3)-f*p(3);

   if y<0.25, f=0; elseif y>0.75, f=1; else f=0.5; end
   p(2)=p0(2)+y*p0(4)-f*p(4);

   h=axes('Position',p);
   setuser(ah,'ih',h);
end

% -------------------------------------------------------------------- %

