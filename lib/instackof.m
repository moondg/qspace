function i=instackof(f,pflag)
% function i=instackof(f [,'-p'])
%    
%    Check whether current position program has the f
%    in its dbstackof in either file or function name.
%    If '-p' is specified, f is interpreted as regexp.
%    
% Wb,Jul05,19

  if nargin>1
       if ~isequal(pflag,'-p'), error('Wb:ERR','\n   ERR invalid usage'); end
  else f=['^' f '\>'];  % note that file name may have extra `.m'
  end

  S=dbstack; i=0;

  for k=1:numel(S)
      if ~isempty(regexp(S(k).file,f)), i=1; break
      elseif ~isempty(regexp(S(k).name,f)), i=1; break
      end
  end

end

