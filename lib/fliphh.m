function fliphh(varargin)
% function fliphh([ah,] lh [,opts])
%
%   flip order of handles lh given the parent handle ah (default: gca)
%
% Options
%   'idx',..  actual index order (rather than n:-1:1)
%
% Wb,Oct29,10

% tags: line handle order

  if ~nargin
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  if nargin>1 && numel(varargin{1})==1 ...
  && ishandle(varargin{1}) && all(ishandle(varargin{2}))
     ah=varargin{1};
     lh=varargin{2}; varargin(1:2)=[];
  else
     ah=gca;
     lh=varargin{1}; varargin(1)=[];
  end

  getopt('init',varargin);
    idx=getopt('idx',[]);
  getopt('check_error');

  ch=get(ah,'children'); lh=lh(:);
  [ih,ip,I]=matchvec(lh,ch,'-s');

  if ~isempty(I.ix1), wbdie(...
     'not all handles found in children of parent\n'); end

  if isempty(idx), idx=numel(ih):-1:1; end

  ch(ip)=ch(ip(idx));
  set(ah,'Children',ch);

end

% -------------------------------------------------------------------- %

