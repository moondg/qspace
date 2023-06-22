function label(varargin)
% Function: label() - adds x/y label and title to current figure
% Usage: label ([ah,] 'xstr', 'ystr', 'tstr' [, yorientation]);
%
%    If axes handles ah are present, apply labeles to every axis.
%
% Wb,Feb28,01

  if nargin && ~isempty(varargin{1}) && all(isaxis(varargin{1}))
       aflag=1; ah=varargin{1}; varargin(1)=[];
  else aflag=0; ah=gca; end

  nargs=numel(varargin); roty=0;
  if nargs, q=varargin{end};
     if isnumeric(q), roty=varargin{q}; 
     elseif ischar(q) && isequal(lower(q),'roty'), roty=1;
     else q=[]; end
     if ~isempty(q), varargin(end)=[]; nargs=nargs-1; end
  end

  if nargs<1 || roty && nargs<2 || nargs>3, helpthis
     if nargin || nargout, wberr('invalid usage'), end
     return
  end

  lh=cell(1,3);
  tags={ 'XLabel','YLabel','Title' };

  for i=1:nargs
     lh{i}=label_set(ah,tags{i},varargin{i});
  end

  if roty, set(lh{2},'Rotation',0); end

end

% -------------------------------------------------------------------- %
function lh=label_set(ah,tag,str)

  lh=get(ah,tag); if isempty(str) && ~ischar(str), return; end
  if iscell(str), str=sprintf(str{:}); end

  if length(ah)>1, lh=cat(1,lh{:}); end
  set(lh,'String',str); % ,'FontSize',fs

end

% -------------------------------------------------------------------- %
