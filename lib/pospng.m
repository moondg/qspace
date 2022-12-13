function h=pospng(fname,x,varargin)
% Function pospng(fname, pos [,opts])
%
%    position png image given in file fname
%    at positon pos relative to current axes
%
% Options
%
%   '-transp'  make white background transparent
%   '-rm'      remove all teximage objects first.
%   'sfac',..  factor on size
%   'tex',..   use output of tex into fname as png.
%
% Examples
%
%    pospng('',[.03,.8],'-rm','sfac',[.95 1.05],'tex',[...
%      'T_K \equiv \sqrt{U\Gamma/2} ' ...
%      'e^{\pi\frac{\epsilon_d(\epsilon_d+U)}{2U\Gamma}}'])
%
% see also teximage.m
% Wb,Aug05,08

  tag='teximage';

  getopt('init',varargin);
     vflag =getopt('-v');
     transp=getopt('-transp');
     sfac  =getopt('sfac',1);
     texs  =getopt('tex','');
     if getopt('-rm'), delete(findall(gcf,'tag',tag)); end
  getopt('check_error');

  if nargin<2
     eval(['help ' mfilename]);
     if nargin || nargout, wberr('invalid usage'), end, return
  end

  if ~isempty(texs)
     if isempty(fname), fname=tempname; end
     if vflag,o={'-v'}; else o={}; end
     latex2png(texs,fname,o{:}); % '-v'
  end

  if isempty(findstr(fname,'.png')), fname=[fname '.png']; end

  ah=gca; fh=gcf; u=get(ah,'Units');
  set(ah,'Units','pixels'); p=get(ah,'Position');
  set(ah,'Units',u);

  I=imread(fname); I=I(end:-1:1,:,:); s=size(I);
  if ~isempty(texs), delete(fname); end

  if numel(sfac)==1, sfac=repmat(sfac,1,2); end
  p = [p(1)+p(3)*x(1), p(2)+p(4)*x(2), sfac.*[s(2)*0.6, s(1)*0.7]];

  ax=axes('visible','off','units','pixels','position',p,'tag',tag);
  h=image('cdata',I,'parent',ax,'tag',tag);

  try
    if transp
      a=I(:,:,3); a=1-a/max(a(:));
      set(h,'alphadata',a,'AlphaDataMapping','scaled');
      if ~isequal(get(fh,'renderer'),'opengl')
	     set(fh,'renderer','opengl');
      end
    else
      set(h,'alphadata',1);
    end
  end

  set(gcf,'CurrentAxes',ah);

end

