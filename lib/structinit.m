function [s,s1]=structinit(varargin)
% function [s,s1]=structinit(varargin)
%
%    init empty structure with given fields (1x0 struct).
%    s1 is an instantiation (1x1 struct) if required.
%
% Wb,May13,09

  m=1; n=0; narg=length(varargin);

  if narg && isnumeric(varargin{1}) && isscalar(varargin{1})
  n=varargin{1}; varargin(1)=[]; narg=narg-1; end

  if narg && isnumeric(varargin{1}) && isscalar(varargin{1}), m=n;
  n=varargin{1}; varargin(1)=[]; narg=narg-1; end

  if narg==1 && isstruct(varargin{1})
     s=reshape(fieldnames(varargin{1}),1,[]); s(2,:)={[]};
     s=repmat(struct(s{:}),m,n);
     return
  end

  e=0; for i=1:narg, if ~ischar(varargin{i}), e=1; break; end, end
  if e
     eval(['help ' mfilename]);
     if nargin || nargout, wberr('invalid usage'), end, return
  end

  if narg
     ss=varargin; ss(2,:)={ cell(m,n) };
     s=struct(ss{:});
     if nargout>1, ss(2,:)={ cell(1,1) }; s1=struct(ss{:}); end
  else
     s=repmat(struct,m,n);
     if nargout>1, s1=repmat(struct,1,1); end
  end

end

