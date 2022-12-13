function mvlabel(varargin)
% function mvlabel([ah,] lb, dx)
%
%    movel label in figure axis, where
%    ah  = axis handle
%    lb  = 'X|Y|Z|T':  x label, y label, z label, title respectively
%    dx = shift in normalized units
%
% Wb,May22,03

  if nargin && all(ishandle(varargin{1}))
     if ~any(isaxis(varargin{1}))
        mv_text_label(varargin{:});
        return
     end
 end

  getopt('INIT',varargin);
     pos =getopt('pos',[]);
     POS =getopt('POS',[]);
  varargin=getopt('get_remaining');

  if ~isempty(POS)
     if ~ischar(varargin{1}), wberr(...
        'invalid usage (no axis handle with POS)'); end
     ah=findall(gcf,'type','axes'); n=numel(ah);
     for i=1:n
        if ~isempty(regexp(get(ah(i),'tag'),'legend|scribe')), continue; end
        mvlabel(ah(i),'pos',POS,varargin{:});
     end
     return
  end

  narg=length(varargin);

  if isempty(pos)
     if narg==2, [ah,lb,dx]=deal(gca,varargin{1:2});
     else
        if narg~=3, eval(['help ' mfilename]); return; end
        [ah,lb,dx]=deal(varargin{1:3});
     end
  else
     if narg==1, ah=gca; lb=varargin{1};
     elseif narg==2, [ah,lb]=deal(varargin{1:2});
     else eval(['help ' mfilename]); return; end
     dx=pos;
  end

  ne=length(dx); unorm=1;
  if ne<1 || ne>3, eval(['help ' mfilename]); return; end

  switch lower(lb)
    case 'x',  h=get(ah,'xlabel'); if ne==1, dx=[0 dx 0]; end
    case 'y',  h=get(ah,'ylabel'); if ne==1, dx=[dx 0 0]; end
    case 'z',  h=get(ah,'zlabel'); if ne==1, dx=[dx 0 0]; end
    case 't',  h=get(ah,'title' ); if ne==1, dx=[0 dx 0]; end
  otherwise
    unorm=0;
    if isnumeric(lb) && isscalar(lb) && ...
       ishandle(lb) && isequal(get(lb,'type'),'text')

       h=lb; dx=[dx,0];
    else
       helpthis, if nargin || nargout, wberr('invalid usage'), end
       return
    end
  end

  if iscell(h)
  h=cat(1,h{:}); end

  if length(dx)==2, dx=[dx 0]; end

  if ~isempty(pos)
      set(h,'Units','norm','Pos',dx);
  else
      for i=1:numel(h)
         u=get(h(i),'Units');
         set(h(i),'Units','norm');
         set(h(i),'Pos',get(h(i),'Position')+dx);
         if ~unorm, set(h(i),'Units',u); end
      end
  end

end

% -------------------------------------------------------------------- %

function mv_text_label(th,dx,varargin)
   set(th,'Units','normalized');
   for i=1:numel(th)
      p=get(th(i),'Position'); p(1:2)=p(1:2)+dx;
      set(th(i),'Position',p);
   end
end

% -------------------------------------------------------------------- %

