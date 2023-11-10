function [dd,qq]=getzdim(A,varargin)
% function [dd,qq]=getzdim(A [,dim,opts])
%
%   get multiplet dimensions of QSpace A.
%
% Options
%
%   -p   product of multiplet dimensions
%   -x   expand to match data (repmat)
%   -u   get dimension of unique qlabel set (returned as 2nd output argument)
%
%  The specification of the dimension (dim) may also take the following
%  string values:
%
% '-op'  gets dimension of IROP set
%  'op'  gets dimension of state space dimension for scalar operator
%        (together with '-u' this is equivalent to 'op') // deprecated
%
% Wb,Jan24,11

  if nargin<1 || numel(A)~=1
     eval(['help ' mfilename]);
     if nargin || nargout, wbdie('invalid usage'), end, return
  end

  getopt('init',varargin);
     pflag=getopt('-p');
     xflag=getopt('-x');
     uflag=getopt('-u');
  dim=getopt('get_last',[]);

  if ischar(dim)
     if isequal(dim,'-op'), uflag=1; end
  end

  if ~uflag
     dd=getzdim_each(A,dim,pflag,xflag);
  else
     if xflag, wblog('WRN','option ''-x'' will be ignored'); end
     [dd,qq]=getzdim_unique(A,dim,pflag);
  end

end

% -------------------------------------------------------------------- %
function dd=getzdim_each(A,dim,pflag,xflag)

  if isempty(A.Q), dd=[]; return; end
  if isempty(A.info)
     dd=repmat(1,[size(A.Q{1}), numel(A.Q)]);
     ns=size(A.Q{1},2);
  else
     r=numel(A.Q); n=numel(A.data); ns=size(A.info.cgr,2);

     dd=zeros(r,n,ns);

     for i=1:n
     for j=1:ns
         s=cgr_size(A,i,j); s(end+1:r)=1; dd(:,i,j)=s;
     end
     end

     dd=permute(dd,[2 3 1]);
  end

  if isempty(dim)
     if pflag || xflag
        wbdie('invalid usage (ignoring -p and -x without dim)'); end
     return
  end

  if ischar(dim)
     if isequal(dim,'op')
        if r~=2 || norm(diff(dd,[],3))>1E-12, wblog('WRN',[...
          'expecting scalar operator with ''op'' !??\n' ... 
          'hint: did you mean ''-op''?']); end
        dim=1;
     else dim, wbdie('invalid usage'); end
  end

  dd=dd(:,:,dim);

  if pflag, dd=prod(dd,2); end
  if xflag
     dd=mat2cell(dd,ones(1,size(dd,1)),size(dd,2));
     for i=1:numel(A.data)
        dd{i}=repmat(dd{i},size(A.data{i},dim),1);
     end
     dd=cat(1,dd{:});
  end

end

% -------------------------------------------------------------------- %
function [dd,qq]=getzdim_unique(A,dim,pflag,xflag)

  if isempty(A.Q), dd=[]; qq=[]; return; end

  r=numel(A.Q);

  if isempty(A.info)
       ns=size(A.Q{1},2);
  else ns=size(A.info.cgr,2); end

  if isempty(dim), idx=1:r;
  elseif isequal(dim,'-op') || isequal(dim,'op')
     if r==2
        dd=ones(1,ns); qq=zeros(1,ns);
        return
     end
     if r~=3, wbdie('invalid operator'); end
     if norm(diff(A.Q{3},[],1))>1E-12
        wblog('WRN','got non-irop (i.e. reducible operator) !??'); end
     idx=3;
  else idx=dim;
  end

  r2=numel(idx); n=size(A.Q{1},1);

  [qq,ii]=uniquerows(cat(1,A.Q{idx}),'-1');
  I=repmat((1:n)',1,r2); I=I(ii);
  J=repmat(1:r2,n,1);    J=J(ii);

  n=numel(ii);

  if isempty(A.info)
     dd=ones(size(qq,1),ns);
  else
     dd=zeros(n,ns);
     for i=1:n
        for j=1:ns
            s=cgr_size(A,I(i),j); s(end+1:r)=1;
            dd(i,j)=s(idx(J(i)));
        end
     end
  end

  if pflag, dd=prod(dd,2); end

end

% -------------------------------------------------------------------- %

