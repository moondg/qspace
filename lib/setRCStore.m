function d=setRCStore(varargin)
% function d=setRCStore([paths])
%
%    Set or inspect RC_STORE environmental variable.
%
% Wb,Aug29,16

  if ~nargin, d=getenv('RC_STORE');
     if isempty(d), d={};
     elseif ~isempty(regexp(d,':')), d=strread(d,'%s','whitespace',':');
     else d={d}; end
     if ~nargout
        print_path('RC_STORE',d);
        clear d
     end
     return
  end

  p0=getenv('RC_STORE');
  p2=strread(strjoin(varargin,':'),'%s','whitespace',':');
  if numel(unique(p2))~=numel(p2)
     wbdie('invalid usage (input paths not unique)');
  end
  for i=1:numel(p2)
     if ~isdir(p2{i}), wbdie('invalid path %s',p2{i}); end
  end

  pstr=strjoin(p2,':');
  if ~isequal(pstr,p0)
     print_path('RC_STORE',p2);
     setenv('RC_STORE',pstr);
  end

end

% -------------------------------------------------------------------- %

function print_path(tag,pp)
  if isempty(pp)
     fprintf(1,'\n   (%s not set)\n\n',tag);
  elseif numel(pp)==1
     fprintf(1,'\n   %s = %s\n\n',tag,pp{:});
  else
     fprintf(1,'\n   %s = \n',tag);
     fprintf(1,'     %s\n',pp{:});
     fprintf(1,'\n');
  end
end

% -------------------------------------------------------------------- %

