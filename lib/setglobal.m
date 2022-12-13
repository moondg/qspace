function setglobal(varargin)
% Function setglobal('var1','var2',...)
%
%   set input arguments global in current workspace
%
% Wb,Feb02,05

% NB! isglobal is an obsolete function and will be
% discontinued (see help isglobal) use this instead:
% if ~isempty(whos('global','variable')), ... 

% NB! iff local and global instance coexist
% clear global does not clear local instance

% fflag = force local instance if it exists global without warning
  fflag=0; if isequal(varargin{1},'-f')
  fflag=1; varargin(1)=[]; end

  if ~isempty(whos('global','sgFlags'))
  wblog('WRN','sgFlags already exists'); end;

  evalin('caller','global sgFlags'); global sgFlags

  nv=length(varargin);
  for iv=1:nv, vn=varargin{iv};

    evalin('caller', sprintf(...
    'sgFlags.il=whos(''%s''); sgFlags.ig=whos(''global'',''%s'');', vn,vn));

    if isempty(sgFlags.il)
    evalin('caller',['global ' vn]); continue;
    elseif sgFlags.il.global, continue; end

    if ~isempty(sgFlags.ig) && ~fflag
       fprintf(1,['\nWarning: %s - variable %s coexists ' ...
       '[local] and global!!\nWarning: %s\n\n'],mfilename,vn,lineno('all'));
    end

    evalin('caller', sprintf(...
    'var__=%s; clear %s; global %s; %s=var__; clear var__',...
    vn,vn,vn,vn));

  end

  clear global sgFlags

end

