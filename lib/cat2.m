function dd=cat2(ic,varargin)
% function dd=cat2(d [,opts,{value}])
%
%    Similar to cat, but adjust dimensions if required
%
%    last argument if of cell type {value} sets value
%    for space added  (default: 0)
%
% Wb,Jul23,09

  if ~nargin || numel(ic)~=1 || any(round(ic)~=ic) || ic<1 || ic>9
     helpthis, if nargin || nargout, wberr('invalid usage'), end
     return
  end

  x=0;
  if iscell(varargin{1})
     x={[]};
  elseif iscell(varargin{end}) && ~iscell(varargin{1})
     q=varargin{end}; 
     if numel(q)==1 && isnumeric(q{1}), x=q{1}; varargin(end)=[]; end
  end

  n=numel(varargin); if ~n, dd=[]; return; end
  s=cell(1,n); rr=zeros(1,n);
  for i=1:n
     if isempty(varargin{i}), varargin{i}=x; end
     s{i}=size(varargin{i}); rr(i)=length(s{i});
  end

  r=max(ic,max(rr)); ss=ones(n,r);
  for i=1:n, ss(i,1:rr(i))=s{i}; end

  smin=min(ss,[],1); smin(ic)=-1;
  smax=max(ss,[],1); smax(ic)=-1;

  if ~isequal(smin,smax)
     for i=1:n, smax(ic)=ss(i,ic);
        if any(ss(i,:)<smax), q=repmat(x,smax);
           if all(ss(i,:))
              s=cell(1,r); for j=1:r, s{j}=1:ss(i,j); end
              s=struct('type','()','subs',{s});
              varargin{i}=subsasgn(q, s, varargin{i});
           else
              varargin{i}=q;
           end
        end
     end
  end

  dd=cat(ic,varargin{:});

end

