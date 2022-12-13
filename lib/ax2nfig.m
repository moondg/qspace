function ah=ax2nfig(varargin)
% function ah=ax2nfig([ah0,][opts])
%
%    takes axis handle given (ah0; or current axis if not specified) and
%    copies the whole graph into a new figure;
%
% Options
%
%   -l   also copy legend axis set
%   -h   also copy header axis set,
%        while reducing the size of the actual axis set to be copied
%
% Wb,Apr03,01 ; Wb,Aug14,11

  if nargin && isaxis(varargin{1})
       ah0=varargin{1}; varargin(1)=[];
  else ah0=gca; end

  getopt('init',varargin);
     hflag=getopt('-h');
     lflag=getopt('-l');
     cmap =getopt('-cm');
     zoom =getopt('zoom',0);
  getopt('check_error');

  set(gca,'selected','off')
  saveDisp(ah0);

  pos=get(groot,'DefaultAxesPosition');
  if hflag && ~zoom, zoom=0.9; end
  if zoom
     pos=[ pos(1:2)+((1-zoom)/2)*pos(3:4), zoom*pos(3:4) ];
  end

  f=figure;
  f0=get(ah0,'parent');

  if hflag
     h=findall(f0,'tag','frame');
     if ~isempty(h)
        h2=copyobj(h,f);
     end
  end

  ah=copyobj(ah0,f);
  set(ah,'Units','Normalized','Position',pos);

  restoreDisp(ah0);
  restoreDisp(ah);

  setax(ah);

  if cmap
     colormap(get(f0,'ColorMap'));
     clim(get(ah0,'CLim'));
  end

  if lflag
    lh=findleg(ah0);
    if ~isempty(lh)
       lp=get(lh,'Position');
       p0=get(ah0,'Position'); l2=ax_pos2pos(lp,p0,pos);
       repax(axes('Position',l2), lh);
    end
  end

  if nargout==0, clear ah; end

end

% -------------------------------------------------------------------- %

function saveDisp(ah)
  for h=findall(ah,'type','line')'
     if isfield(get(h),'DisplayName')
        setuser(h,'disp',get(h,'DisplayName'));
     end
  end
end

function restoreDisp(ah)
  for h=findall(ah,'type','line')'
     if isfield(get(h),'DisplayName')
        set(h,'DisplayName',getuser(h,'disp','-rm'))
     end
  end
end

% -------------------------------------------------------------------- %

