function A=itagrep(A,varargin)
% function A=itagrep(A [,idx], <arguments to regepxrep>)
%
%    Apply given regexprep() pattern to specified itags
%    idx specified a subset of itags to consider;
%    (otherwise, by default, all itags are considered)
%
% Wb,Jan14,15

  if nargin<3 || numel(A)~=1
     helpthis, if nargin || nargout, wberr('invalid usage'), end
     return
  end

  t=getitags(A);

  if isnumeric(varargin{1})
       it=varargin{1}; varargin(1)=[];
  else it=1:numel(t); end

  t(it)=regexprep(t(it),varargin{:});
  A.info.itags=t;

  if ~nargout, vname=inputname(1);
     if ~isempty(vname)
        assignin('caller',vname,A);
        clear A
     end
  end

end

