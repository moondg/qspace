function cgs=cgs2double(varargin)
% function cgs=cgs2double(cgs,...)
% Wb,Aug11,14

  getopt('init',varargin);
     if getopt('-sp'); spflag=1;
     elseif getopt('-x'); spflag=0; else spflag=-1; end
  varargin=getopt('get_remaining'); narg=length(varargin);

  if ~narg || iscell(varargin{1}) && narg>1
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  if iscell(varargin{1})
     if narg==1, varargin=varargin{1};
     else error('Wb:ERR','\n   ERR invalid cell-usage'); end
  end

  n=numel(varargin); cgs=cell(1,n);
  for i=1:n
     cgs{i}=cgs2double_1(varargin{i},spflag);
  end

  if n==1, cgs=cgs{1}; end

end

% -------------------------------------------------------------------- %
 
function cgs=cgs2double_1(cgs,spflag)

  if ~isfield(cgs,'S') || ~isfield(cgs,'data') 
     error('Wb:ERR','\n   ERR invalid MPFR structure');
  end

  if isfield(cgs,'idx')
     for i=1:numel(cgs);
        q=cgs(i).data; if size(q,1)==1, q=q'; end

        if isempty(q), 
           cgs(i).data=zeros(cgs(i).S); continue;
        end

        q=mpfr2dec(q);

        if isempty(cgs(i).S) && isempty(cgs(i).idx)
           s=numel(q); i1=1:s;
           cgs(i).data=sparse(i1,i1,q,s,s); s=s([1 1]);
        else
           s=double(cgs(i).S); cs=[1, cumprod(s)];
           idx=double(cgs(i).idx);
           i1 = idx(:,1:end-1) * cs(1:end-2)'+1;
           cgs(i).data = sparse(i1, idx(:,end)+1, q, cs(end-1),s(end));
        end

        if ~spflag || spflag<0 && length(cgs(i).S)~=2
           cgs(i).data = reshape(full(cgs(i).data), s);
        end
     end
  else
     for i=1:numel(cgs);
        q=cgs(i).data; if size(q,1)==1, q=q'; end
        cgs(i).data= reshape(mpfr2dec(q), cgs(i).S);
     end
  end

  if numel(cgs)==1
       cgs=cgs.data;
  else cgs=reshape({cgs.data},size(cgs));
  end

end

% -------------------------------------------------------------------- %
 
