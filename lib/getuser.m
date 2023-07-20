function v=getuser(varargin)
% Function u=getuser([hh,] tag)
%
%    Get field tag of userdata of current axes (or hh).
%
% Options
%
%    '-rm'  afer obtaining given data, remove field from user data.
%
% Wb,Sep02,08

  if nargin && all(ishandle(varargin{1}))
       ah=varargin{1}; varargin(1)=[];
  else ah=gca;
  end

  rflag=0;
  for i=1:numel(varargin)
     if isequal(varargin{i},'-rm'), rflag=1; varargin(i)=[]; break; end
  end
  narg=numel(varargin);

  nd=numel(ah); 

  if ~narg
     v=get(ah,'UserD');
     if nd>1, v=catc(2,v{:}); end
     return
  end

  if narg~=1 || ~ischar(varargin{1})
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  f=varargin{1}; v=cell(1,nd);

  for i=1:nd, h=ah(i);
     u=get(h,'UserData');
     if ~isempty(h) && isfield(u,f)
        v{i}=getfield(u,f);
        if rflag, set(h,'UserD',rmfield(u,f)); end
     else v{i}=[]; end
  end

  if nd==1, v=v{1};
  elseif nd>1
     s1=size(v{1});
     for i=2:nd
        if ~isequal(size(v{i}),s1), i=0; break; end
     end
     if i, v=[v{:}]; end
  end

end

