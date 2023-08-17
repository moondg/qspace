function hflag=helpthis(nargs_out,varargin)
% function q=helpthis([nargs_out, varargin{:}])
%
%    Show help for the caller function.
%    If this is called via helpthis(nargout, varargin{:})
%    then the help is only shown for -? or -h usage,
%    in which case hflag=1 is returned (hflag=0 otherwise).
%
% Wb,Feb23,11

  if nargin, hflag=true;
     if nargs_out || nargin>2, hflag=false;
     elseif nargin==2 || ~ischar(varargin{1}), hflag=false;
     elseif isempty(regexp(varargin{1},'^(-[h?]|--help)$')), hflag=false;
     end
     if ~hflag, return; end
  end

  I=dbstack(1,'-completenames');
  n=numel(I);
  s={};

  for i=1:n
     if i>1 && ~isequal(I(i).file,I(i-1).file), i=i-1; break; end
     s={I(i).name, s{:}};
  end

  f=I(i).file;
  i=regexp(f,'Class/@');
  if ~isempty(i)
     f=f(i+7:end); % excluding '@' // Wb,Aug23,16
  else
     i=find(f=='/');
     if ~isempty(i), f=f(i(end)+1:end); end
  end
  if isequal(f(end-1:end),'.m'), f=f(1:end-2); end

  n=numel(s); fprintf(1,'\n');
  if n==1, help(f);
  elseif n==2 && ~isequal(f,s{end})
     n=[f '>' s{end}]; help(n);
  elseif n>1, help(f);
     if n>2, fprintf(1,'  NB! got recursive call at level=%g\n',n); end
  else wbdie('invalid usage (don''t know what to do with n=%g)',n); end

end

