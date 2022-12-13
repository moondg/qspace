function B=transform(varargin)
% function B=transform(A,[opts,]U)       % keep U last
%
%    B=U'*A*U
%
% function transform(A,B,...,[opts,]U)   % keep U last
%
%    if more than two input args are specified, 
%    variables are set back in caller space U'*A*U.
%
% Wb,Nov17,11

  getopt('init',varargin);
     sU=getopt('-n',inputname(nargin));
  varargin=getopt('get_remaining'); narg=length(varargin);

  n=narg-1;
  U=varargin{narg}; if isempty(sU), sU='(U)'; end

  save2caller=(~nargout && narg>2);

  for k=1:n
     B=varargin{k};
     for i=1:numel(B)
        if ~isempty(B(i).op)
           s=get_mstr(B(i)); if ~isempty(regexp(s,'[+-]')), s=['[ ' s ' ]']; end
           B(i).istr=[ sU '''*' s '*' sU ];
           B(i).op=U'*B(i).op*U; if ~isempty(B(i).hc)
           B(i).hc=U'*B(i).hc*U; end
        end
     end
     if save2caller
        s=inputname(k); if isempty(s), error('Wb:ERR',...
         '\n   ERR failed to obtain variable name in caller spce'); end
        assignin('caller',s,B);
     end
  end

end

