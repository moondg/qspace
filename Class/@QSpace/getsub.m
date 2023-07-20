function [A,B]=getsub(A,varargin)
% Usage #1: function [A,B]=getsub(A,idx)
%
%    Select given subspace where idx represents direct index
%    into records of QSpace.
%
% Usage #2: function [A,B]=getsub(A,q,dim)
%           function [A,B]=getsub(A,{q} [,dim]) // Wb,Feb21,21
%
%    Select subspace that has given symmetry sector(s) in given dimension(s).
%    if q is a cell array, then this directly relates to to dimensions in A
%    e.g. getsub(ops,{[],[],0}) = getsub(ops,0,3) picks scalar sector
%    on 3rd dimension; here '0' or 0 can be used as a shortcut to zeros(1,m)
%    in q, with m the number of symmetry labels.
%
% Usage #3: function [A,B]=getsub(A,'-0')       // Wb,Dec04,20
%           function [A,B]=getsub(A,'-0!')
%
%    Assuming A is an irop (rank-3 with irop-index on 3rd leg)
%    select scalar operator part (q-labels 0); this is equivalent
%    to usage #2, having getsub(A,[0 0 ..],3));
%
%    in case of '-0!', also remove singleton irop dimension
%    (in case it is not a singleton, i.e., by having multiple
%    q=0 multiplets, an error is issued)
%
% Options
%
%   -d   return A.data (rather than A itself)
%
%   If second output argument is requested, this will contain
%   whatever was left in A_in as compared to what was selected
%   in A_out.
%
% Wb,Aug20,08 ; Wb,Jun19,13

  if nargin<2 || nargin>3
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  sflag=0;
  for i=2:nargin
     if ischar(varargin{i-1}), s=varargin{i-1};
        if regexp(s,'^-0!?')
           if s(end)=='!', sflag=2; else sflag=1; end
           varargin(i-1)=[]; break;
        end
     end
  end

  getopt('init',varargin);
     dflag=getopt('-d');
  args=getopt('get_remaining');

  nA=numel(A);
  if nargout<=1
     for k=1:nA, [A(k)     ]=getsub_1(A(k),sflag,args{:}); end
  else
     B=QSpace(size(A));
     for k=1:nA, [A(k),B(k)]=getsub_1(A(k),sflag,args{:}); end
  end

  if dflag, sA=size(A);
     A=reshape({A.data},sA);
     ss=zeros(1,nA); for i=1:nA, ss(i)=numel(A{i}); end
     if all(ss<=1)
        i=find(ss==0); A(i)={[]};
        A=reshape([A{:}],sA);
     end
  end

end

% -------------------------------------------------------------------- %
function [A,B]=getsub_1(A,sflag,varargin)

  B=QSpace; if isempty(A.Q), return; end
  m=size(A.Q{1},2);

  if sflag
     rA=numel(A.Q); if ~rA || rA==2, return; end
     if rA~=3, wbdie('invalid usage #3 (got rank-%g QSpace)',rA); end

     idx=matchIndex(A.Q{3},zeros(1,m));

  elseif nargin==3 && ~iscell(varargin{1}), idx=varargin{1};
     if ~isnumeric(idx), wbdie('invalid idx'); end
     if isempty(idx)
        if nargout>1, B=A; end
        A=QSpace; return
     end
     if ~isvector(idx)
        wblog('WRN','got non-vector as index (%s)',sizestr(idx));
        idx=idx(:);
     end
     if numel(unique(idx))~=numel(idx)
        wbdie('invalid usage (non-unique idx)');
     end

  else
     if nargin>3, k=varargin{2}; else k=[]; end
     q=varargin{1}; n=numel(q);
     if iscell(q)
        if isempty(k)
           k=1:n; for i=1:n, if isempty(q{i}), k(i)=0; end; end
           k=k(find(k));
        else
           s=[ n, numel(k) ]; if diff(s)
           wbdie('invalid usage #3 (size mismatch %g/%g)',s); end
        end

        qq=[];
        for i=1:n, if isempty(q{i}), continue; end
           if isequal(q{i},'0') || isequal(q{i},0), q{i}=zeros(1,m);
           elseif size(q{i},2)~=m, wbdie(...
             'invalid usage #2 [size(q{%g},2) = %g/%g]',i,size(q{i},2),m);
           end
           if ~isempty(qq)
              i2=vkron(1:size(qq,1), 1:size(q{i},1));
              qq=[ qq(i2(:,1),:), q{i}(i2(:,2),:) ];
           else qq=q{i}; end
        end
        q=qq;

     elseif isempty(k)
        wbdie('invalid usage #2 (missing argument dim)');
     end

     ia=matchIndex(cat(2,A.Q{k}),q);
     idx=unique(ia);
  end

  if nargout>1, B=A;
     if ~isempty(B.Q) && ~isempty(B.Q{1})
        n=size(B.Q{1},1); j=1:n; j(idx)=[];

        if isempty(j), B=QSpace;
        else
           for p=1:length(B.Q), B.Q{p}=B.Q{p}(j,:); end
           B.data=B.data(j);
           if isfield(B.info,'cgr') && ~isempty(B.info.cgr)
              B.info.cgr=B.info.cgr(j,:);
           end
        end
     end
  end

  for p=1:length(A.Q), A.Q{p}=A.Q{p}(idx,:); end
  A.data=A.data(idx);
  if isfield(A.info,'cgr') && ~isempty(A.info.cgr)
     A.info.cgr=A.info.cgr(idx,:);
  end

  if sflag>1
     d=getDimQS(A); if d(end)~=1, wbdie(...
       'invalid usage (got dimension d=%g on scalar sector)',d(end)); end
     A=squeeze(A,3);
  end

end

% -------------------------------------------------------------------- %

