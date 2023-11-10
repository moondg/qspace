function EE=reduceEE(EE,varargin)
% Function EE=reduceEE(EE)
% 
%    remove degenerate levels from nrg_plot
%    also remove levels outside yrange if specified.
% 
% Options
% 
%   'yl',...  remove levels outside given yrange (Inf)
%   'eps',..  remove levels whos absolute difference to its neighboring
%             record is less than (1E-2)
% 
% Wb,Jun09,08

  getopt('init',varargin);
     yl =getopt('yl',Inf);
     eps=getopt('eps',1E-2);
  getopt('check_error');

  if nargin<1 || ~isnumeric(EE)
     eval(['help ' mfilename]);
     if nargin || nargout, wbdie('invalid usage'), end, return
  end

  s=size(EE); if s(1)<s(2), wblog('ERR',...
  'expecting energy spectrum along columns !?? (%s)',vec2str(s)); end

  n=s(1); mark=zeros(n,1); y0=EE(1,:);

  for i=2:n
      yy=EE(i,:); if all(yy>yl), mark(i)=1; continue; end

      j=find(~isnan(yy) & ~isnan(y0) & (y0~=0 | yy~=0));
      if norm(yy(j)-y0(j))<eps, mark(i)=1; end

      y0=yy;
  end

  i=find(mark); m=length(i);
  wblog('<i>',...
     'removing %g degenerate records (%g, E<%g: %.3g%%)',...
      m,eps,yl,100*m/n ...
  )

  EE(i,:)=[];

end

