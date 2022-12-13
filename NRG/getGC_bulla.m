function [ac,I]=getGC_bulla(varargin)
% Function: [ac,I]=getGC_bulla(om, gi,fi [,opts])
% Function: [ac,I]=getGC_bulla([om gi], [om fi] [,opts])
%
%    calculate improved spectral function using Bulla'98
%    with gi and fi preferentially already but not necessarily
%    smoothened (note that KramersKronig within discrete data
%    is also possible, but when smoothening this deteriorates
%    somewhat the sum-rule for spectral function).
%
%    (SIAM) parameters are read from global/base `param'.
%
% Options
%
%   'omfac',..   linear interpolation for energies < min(om)*omfac (2)
%   '-iB'        indicate that Hamiltonian used -B*Sz (note the sign)
%   '--Uh'       shift epsd by U/2 (e.g. when the interaction
%                was written as U*nup*ndown => U*(n-1)^2
%
%    ...         model specific parameters; whatever is left then
%                is forwarded to KKreal.
%
% Wb,Oct22,05 ; Wb,Jan23,11 ; Wb,Dec11,20

% outsourced from $VER//xGF
% see also getGC_JJH.m

  for n=1:nargin
     if ~isnumeric(varargin{n}), n=n-1; break; end
  end

  if ~nargin || n<2 || n>3
     helpthis, if nargin || nargout, wberr('invalid usage'), end
     return
  end

  wusage=1;

  if n==2
     [gi,fi] = deal(varargin{1:2}); varargin(1:2)=[];
     s=[size(gi,2), size(fi,2)];
     if s(1)>1 && ~diff(s)
        om=gi(:,1);
           if length(om)~=size(fi,1) ||  norm(om-fi(:,1))>1E-15
           wberr('gi and fi must have the same omega range'); end
        gi=gi(:,2:end); fi=fi(:,2:end);
        wusage=2;

     elseif s(1)==1 && s(2)>1 && ~mod(s(2),2), m=s(2)/2;
        om=gi; gi=fi(:,1:m); fi=fi(:,m+1:end);
        wusage=2.2;

     else wberr('invalid usage (s=[%g,%g])',s);
     end

  elseif n==3,
     [om,gi,fi] = deal(varargin{1:3}); varargin(1:3)=[];

     if ~isvector(om) || length(om)~=size(gi,1) || length(om)~=size(fi,1)
        wberr('invalid omega (size mismatch)'); end
     wusage=3;

  else wberr('invalid usage'); end

  if ~isvector(om)
     eval(['help ' mfilename]);
     if nargin || nargout, wberr('invalid usage'), end, return
  else om=om(:); end

  global param
  if isempty(param), getbase param; end
  if isempty(param), wberr('ERR failed to get param(eters)'); end

  getopt('init',varargin);
     omfac=getopt('omfac',2);
     if getopt('-q'), vflag=0; else vflag=1; end
     iBflag=getopt('-iB');

     Uhflag=getopt('--Uh');
     ip=getopt('ip',[]);

  oKK=getopt('get_remaining');

  if ~all(isfield(param,{'Gamma','epsd','U'})), wberr('invalid usage'); end

  I=struct;
     for s={'U','Gamma','epsd','B'}
         I=setfield(I,s{1},getfield(param,s{1})); end
     if ~isfield(I,'B'), I.B=0; end
     if iBflag, I.B=-I.B; end
  structexp(I);

  m=size(fi,2);
  if m>2, wberr('invalid input data (%g)',size(fi,2)); end
  if B && mod(m,2) && vflag
     wblog('WRN','assuming spin up for improved spectral function');
  end

  q=[ numel(U), numel(epsd), numel(Gamma) ];
  if isempty(ip)
     if B, ip=floor((2:m+1)/2);
     else  ip=1:m; end
  end
  I.ip=ip; I=struct('param',I);

  if ~isreal(gi), e=intxy(om,imag(gi));
     if e>1E-12, wblog('WRN','skipping imaginary part @ %.3g',e); end
     gi=real(gi);
  end
  if ~isreal(fi), e=intxy(om,imag(fi));
     if e>1E-12, wblog('WRN','skipping imaginary part @ %.3g',e); end
     fi=real(fi);
  end

  gi=(-pi)*gi; gg=zeros(size(gi));
  fi=(-pi)*fi; ff=zeros(size(fi));

  l=max(max(ip),m);
  U(end+1:l)=U(end);
  epsd(end+1:l)=epsd(end);   if Uhflag, epsd=epsd+U/2; end
  Gamma(end+1:l)=Gamma(end);

  eB=[-B/2, +B/2];
  dd=get_delta(om);

  for k=1:m, l=ip(k);
     gg(:,k) = complex( KKreal(om,gi(:,k),oKK{:}), gi(:,k) );
     ff(:,k) = complex( KKreal(om,fi(:,k),oKK{:}), fi(:,k) );

     if U(l)
          Sf=ff(:,k)./gg(:,k);
     else Sf=0; end

     gc(:,k) = 1./( om - epsd(l) - eB(mod(k-1,2)+1) - Gamma(l)*dd - U(l)*Sf);

     US(:,k) = Sf;
  end

  add2struct(I,gg,ff,gc,US,Uhflag,iBflag,wusage);

  if omfac>0
     aom=abs(om); x0=min(aom)*omfac;

        i0=find(aom<x0);
        i1=find(om==max(om(find(om>=-x0))),1);
        i2=find(om==min(om(find(om>=+x0))),1);

        add2struct(I,omfac,x0);

     if isempty(i0) || isempty(i1) || isempty(i2), return; end
     om0=om(i0); x12=[om(i1),om(i2)];

     for k=1:size(fi,2)
        gc(i0,k) = polyval(polyfit(x12,[gc(i1,k),gc(i2,k)],1),om0);
        US(i0,k) = polyval(polyfit(x12,[US(i1,k),US(i2,k)],1),om0);
     end
  end

  ac=-(1/pi)*imag(gc);

end

% -------------------------------------------------------------------- %

function dd=get_delta(oz)

  if isreal(oz), oz=complex(oz, +1E-12); end

  dd = (1/pi)*log((oz+1)./(oz-1));

end

% -------------------------------------------------------------------- %

