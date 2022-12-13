function S=rmfields(S,varargin)
% Function: rmfields(S, gpattern1, ...)
%
%    remove fields with given grep pattern
%
% Wb,Dec16,07

  if nargin<2 || ~isstruct(S)
     helpthis, if nargin || nargout, wberr('invalid usage'), end
     return
  end

  f=fieldnames(S); nf=length(f);

  for i=1:numel(varargin), if ~ischar(varargin{i})
     if iscell(varargin{i}), wberr([...
       'invalid usage: got cell array' ... 
       'hint: use standard rmfield() instead.']);
     else wberr('invalid usage: string (g)pattern(s) expected!'); end
  end, end

  gpat=varargin; gpat(2,1:end-1)={'|'}; gpat=cat(2,gpat{:});

  mark=ones(1,nf);
  for i=1:length(f)
  if isempty(regexp(f{i},gpat)), mark(i)=0; end, end

  f=f(find(mark)); if isempty(f), return; end
  S=rmfield(S,f);

return

  n=inputname(1);
  if nargout==0 && ~isempty(n)
     assignin('caller',n,S);
     clear S
  end

end

