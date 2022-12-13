function label(varargin)
% Function: label() - adds x/y label and title to current figure
% Usage: label ([ah,] 'xstr', 'ystr', 'tstr' [, yorientation]);
%
%    If axes handles ah are present, apply labeles to every axis.
%
% Wb,Feb28,01

  if numel(varargin) && ~isempty(varargin{1}) && all(isaxis(varargin{1}))
       aflag=1; ah=varargin{1}; varargin(1)=[];
  else aflag=0; ah=gca; end

  roty=0;
  if numel(varargin)
     if isnumeric(varargin{end})
        roty=varargin{1}; varargin(end)=[];
     elseif isequal(lower(varargin{end}),'roty')
        roty=1; varargin(end)=[];
     end
  end

  narg=numel(varargin);
  if narg<1 || roty && narg<2 || narg>3
     eval(['help ' mfilename]);
     if nargin || nargout, wberr('invalid usage'), end, return
  end

  tags = { 'XLabel','YLabel','Title' }; lh=cell(1,3);

  for i=1:narg
     lh{i}=label_set(ah,tags{i},varargin{i});
  end

  if roty, set(lh{2},'Rotation',0); end

end

% -------------------------------------------------------------------- %
function lh=label_set(ah,tag,str)

  lh=get(ah,tag); if isempty(str) && ~ischar(str), return; end

  if length(ah)>1, lh=cat(1,lh{:}); end
  set(lh,'String',str); % ,'FontSize',fs

end

% -------------------------------------------------------------------- %
