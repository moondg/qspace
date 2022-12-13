function [i,x]=isset(vn,i0)
% function [i,val]=isset(vn [,i0])
%
%    check whether variable with name 'vn' is defined in current workspace.
%    if variable is not found, then check value of i0, instead, if specified.
%
%    if vn is numeric, this returns 1 if vn is non-empty and != 0.
%
% Wb,Feb03,10

% tags: gotq exists

  if nargin==1 && isnumeric(vn)
     if ~isempty(vn) && (numel(vn)>1 || vn), i=1; else i=0; end
     return 
  end

  if nargin<1 || nargin>2
     helpthis, if nargin || nargout, wberr('invalid usage'), end
     return
  end

  if ~iscell(vn), vn={vn}; end; nv=numel(vn);
  if nargin<2, i0=[]; end

  qx=cell(nv,2);
  for i=1:nv, v=vn{i};
     if ~ischar(v) || isempty(v), vn
        wbdie('invalid variable name (strings expected)'); end

     q={'g_isset__','(undef)'};
     setuser(groot,q{:}); evalin('caller',[ 'if exist(''' v ...
        ''',''var''), setuser(groot,''' q{1} ''', ' v '); end']);
     x=getuser(groot,q{1},'-rm');

     if isequal(x,q{2}), x=i0; end
     if isempty(x) || isequal(x,0) || isequal(x,'0') && ischar(x)
     qx{i,1}=0; else qx{i,1}=1; end
     qx{i,2}=x;
  end

  if nv==1, x=qx{2}; i=qx{1}; 
  else
     x=qx(:,2)'; i=[qx{:,1}];
  end

end

