function sms(varargin)
% Function sms([[hh,] msize] [,OPTS])
% set marker size
%
%    hh  axis and/or specific line handle(s)  {current axis}
%        fill markers and set their size {6}
%
%    remaining OPTS be applied to all line handles
%    prior to setting marker size.
%
% SMS options
%
%   '-fix' fix markers for existing line handles
%          eg. if color changed, also change facecolor
%
%    remaining options are directly set to line handles.
%
% Example: sms(4,'Marker','o')
%
% Wb,Apr03,01 Wb,Dec18,07

  getopt('init',varargin);
     dofix=getopt('-fix');
  varargin=getopt('get_remaining'); narg=length(varargin);

  if nargin<1 || narg && ~isnumeric(varargin{1}) && ~all(ishandle(varargin{1}(:)))
     eval(['help ' mfilename]);
     if nargin || nargout, wbdie('invalid usage'); end
     return
  end

  lh=[]; ah=[]; msize=[];
  if narg
     if all(isline(varargin{1}))
        lh=varargin{1}; varargin(1)=[]; narg=numel(varargin);
     elseif all(isaxis(varargin{1}))
        ah=varargin{1}; varargin(1)=[]; narg=numel(varargin);
     end
  end
  if narg && isnumber(varargin{1})
     msize=varargin{1}; varargin(1)=[]; narg=numel(varargin);
  end

  if isempty(ah) && ~dofix, ah=gca; end
  if isempty(lh), lh=findall(ah,'Type','line'); end

  if narg && ~ischar(varargin{1}), {msize,dofix, varargin{:}}
     if nargin || nargout
        wbdie('invalid usage'); else helpthis; end
     return
  end

  if dofix
     if isempty(lh), lh=findall(gca,'Type','line'); end
     for h=reshape(lh,1,[])
        c=get(h,'color');
        set(h,'MarkerEdgeColor',c,'MarkerFaceColor',c)
     end
     if ~narg, return; end
  end

  if ~isempty(msize)
     for h=reshape(lh,1,[]), c=get(h,'color');
        set(h, varargin{:}, ...
       'MarkerSize', msize,'MarkerEdgeColor',c,'MarkerFaceColor',c);
     end
  end

end

