function s=cgr_size(A,i,j)
% function s=cgr_size(A,i,j)
% Wb,Jan08,12

% e.g. used by display()

  if nargin~=3 || numel(A)~=1 || numel(i)~=1 || numel(j)~=1
     error('Wb:ERR','\n   ERR invalid usage'); end

  r=numel(A.Q);
  if ~gotCGS(A), s=ones(1,r); return; end

  c=A.info.cgr;
  if iscell(c), c=c{i,j}; else c=c(i,j); end

  if isfield(c,'size'), s=c.size;
  elseif isnumeric(c),  s=size(c);
  else wberr('invalid cgr data'); end

  l=numel(s);
  if l<r, s(end+1:r)=1;
  elseif l==r+1
     if isfield(c,'cgw'), w=size(c.cgw);
        if numel(w)~=2 || diff(w)>0 || any(w>s(end)), wberr(...
            'unexpected OM dimensions of cgr(%d,%d).cgw @ OM=%d',i,j,s(end));
        else s(end)=min(s(end),w(2));
        end
     end
  elseif l>r+1
     wberr('unexpected cgr rank %d/%d',l,r);
  end

end

