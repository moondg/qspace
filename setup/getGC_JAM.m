function [gc,I,lh]=getGC_JAM(om,G,F,varargin)
% Function GC=getGC_JAM(om,G,F)
%
%   calculate improved spectral function using Bulla'98
%   copy/paste from xGF2 routine
%   adapted for hybrid SIAM with U=epsd=0, but finite (negative) J
%   [see discussion with A. Rosch]
%
% Example: getGC_JAM(ox, ax(:,1), ax(:,2))
%
% Wb,Aug27,07

  global param

  getopt('init',varargin);
     tflag=getopt('-t') | getopt('tst');
     qflag=getopt('-q');
     rmflg=getopt('-rm');
     delta=getopt('delta',[]);
  getopt('check_error');

  if rmflg || nargin && isequal(om,'-rm')
     delete(findall(0,'type','line','tag','GC_JAM'));
     if nargin==1, return; end
  end

  if nargin<3
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

  if isempty(param) || ~isfield(param,'J') || ~isfield(param,'Gamma')
     eval(['help ' mfilename]);
     error('Wb:ERR','field J and Gamma required in param')
  end

  structexp(param)

  if exist('U','var') && U~=0 || exist('epsd','var') && epsd~=0
     eval(['help ' mfilename]);
     warning('Wb:WRN','param: finite U or epsd - will be ignored')
     inl 1; disp(param);
  end

  if ~qflag
  if tflag==0
     h=findall(0,'type','figure','tag','rsmoothSpec');
     if isempty(h), eval(['help ' mfilename]); return, end
     figure(h);

     h1=findall(gcf,'type','axes','tag','Spec01');
     h2=findall(gcf,'type','axes','tag','Spec02');
     if isempty(h1) || isempty(h2), eval(['help ' mfilename]); return, end
     ah=[h1,h2];
  else
     ah=smaxis(1,1,'tag',mfilename);
     setax(ah);
  end
  end

% -------------------------------------------------------------------- %

  if ~isempty(delta)
     G=fixZeroData(om,G,delta);
     F=fixZeroData(om,F,delta);
  end

  gr = KKreal(om,G);
  fr = KKreal(om,F);

  if ~isempty(delta)
     gr=fixZeroData(om,gr,delta);
     fr=fixZeroData(om,fr,delta);
  end

  gg = complex(gr,G);
  ff = complex(fr,F);

  SJ = J*(ff./gg);
  Dl = DeltaG(om,Gamma);
  gc = -(1/pi)./( om - epsd - Dl - SJ);

  afac=pi*Gamma;

  if nargout>1
  I=add2struct('-',gg,ff,SJ,Dl); end

  if ~qflag
  if tflag==0
     delete(findall(gcf,'type','line','tag','GC_JAM'));

     set(gcf,'CurrentAxes',ah(1)); hold on
     h1=plot(om,afac*imag(gc),'r','tag','GC_JAM');
     yl=ylim; if yl(2)<1, yl(2)=1.1; ylim(yl); end

     set(gcf,'CurrentAxes',ah(2)); hold on
     h2=plot(abs(om),afac*imag(gc),'r','tag','GC_JAM');
     yl=ylim; if yl(2)<1, yl(2)=1.1; ylim(yl); end

     lh=[h1,h2];

  else
     plot(om,afac*G, om,afac*F); hold on
     lh=plot(om,afac*imag(gc),'r');
     xlim([-1 1]*0.2)
  end
  end

  if ~nargout, clear gc; end

end

% -------------------------------------------------------------------- %

function dd=DeltaG(oz, Gamma)

  if isreal(oz), oz=complex(oz, +1E-12); end

  dd = (Gamma/pi)*log((oz+1)./(oz-1));

end

% -------------------------------------------------------------------- %

