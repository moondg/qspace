function [A,varargout] = QSpace(varargin)
% function [A,varargout] = QSpace(...)
%
%    Construct QSpace class object.
%
%    This serves to wrap QSpace structures into QSpace objects,
%    and also can be used to generate several specific QSpace objects.
%    For historical reasons, this also includes low-level QSpace
%    contructors for the case of only abelian quantum numbers.
%
% Usage:
%
% #1: A = QSpace(Q,data)
%     A = QSpace(A)   conversion e.g. from struct
%     A = QSpace(A,'check')  also performs full check on QSpace structure
%
% #2: A = QSpace(q1, data1, q2, data2, ...) % see initQSpace
%
% #3: A = initQSpace(Q1,M,DB)
%         square operator M with row/col Q1 block size given by DB
%
% #4: A = initQSpace(Q1,Q2,M,'map')
%         matrix M with row/col QIDX given by Q1/Q2
%         eg. used for projectors
%
% #5: A = QSpace(Q,M,'operator');
% #6: A = QSpace(A,'unity');
%
% #7: A0=initQSpace(QL,Qs1,Qs2,...); setup identity tensor
%         for product spaces (QL,Qs1,...)
%  
%         it must hold: size(QL,2) = size(Qs1,2) = size(Qs2,2) = ... = QDIM.
%         Qs and QL span the full Hilbert space, hence certain quantum
%         lables may appear several times.
%
% Wb,May12,06 ; Wb,Apr16,16

%  for testing purposes
%  8: A = QSpace(CData);  % generate a QSpace based on a particular CData

  if nargin==1, A=varargin{1};
     if isa(A,'QSpace')
        return
     elseif isstruct(A) && isfield(A,'data')
        A=struct_to_QSpace(A); return
     elseif isempty(A)
        A=class(get_struct(),'QSpace'); return
     elseif isequal(A,'--dummy');
        A=class(get_struct({1,1},{1},{'','*'}),'QSpace'); return
     elseif iscell(A)
        na=numel(varargin{1});
        if na && ischar(A{1})
           A=class(get_struct(),'QSpace');
           A.info.itags=varargin{1}; return
        else
           clear A
           for k=na:-1:1, A(k)=varargin{1}{k}; end
           A=reshape(A,size(varargin{1}));
           if ~isa(A,'QSpace')
               A=struct_to_QSpace(A); 
               check_QSpace(A);
           end
           return
        end
     end
  elseif ~nargin, A=class(get_struct(),'QSpace'); return
  end

  if isnumeric(varargin{1})
     q=zeros(nargin,2);
     for i=1:nargin
        q(i,:)=[isnumeric(varargin{i}), numel(varargin{i})];
        if ~q(i), i=-1; break; end
     end
     if i==nargin && (i==1 || all(q(:,2)<=1))
        i=find(~q(:,2)); if ~isempty(i), wbdie(...
          'invalid usage (arg #%d contains empty value for dimension)',i(1));
        end
        d=[varargin{:}];
        A=repmat(class(get_struct(),'QSpace'),d);
        return
     end
  end

  done=0;

  if nargin==1 || nargin==2 && isequal(varargin{2},'check')
     if isa(A,'QSpace'), done=1;
     elseif isstruct(varargin{1}), A=varargin{1};
        if isfield(A,'data')
           A=struct_to_QSpace(varargin{1}); done=1;
        elseif isfield(A,'cgd') && isfield(A,'qset')
           A=initQSpaceCD(A,varargin{2:end}); done=1;
        else wbdie('invalid usage'); end
     end
     if done
        if nargin>1, check_QSpace(varargin{1}); end
        return
     end
  end

  varargout=cell(1,nargout-1);

% usage #1: may specify QSpace with single block entry
% with Q-labels collected into rows for initialization
% e.g. wit NO symmetries e.g. for testing purposes
% and also for consistency with older setups // Wb,Sep14,18
%
% Example: E=QSpace([0;0],eye(4),{'','*'}); % qq, data, info
% where info must either be a cell or struct (QSpace/get_struct.m)

  if nargin==3 && isnumeric(varargin{1}) && ~ischar(varargin{3})
     if size(varargin{1},1)>1
        q=varargin{1}; s=size(q);
        varargin{1}=mat2cell(q,ones(1,s(1)),s(2))';
        if ~iscell(varargin{2}), varargin{2}={varargin{2}}; end
     end
  end

  if iscell(varargin{1}) || isstruct(varargin{1}), A=varargin{1};

     if nargin==2
          A=get_struct(varargin{1}, varargin{2}(:), []);
     else A=get_struct(varargin{1}, varargin{2}(:), varargin{3});
     end

     if ~isempty(A.Q) && size(A.Q{1},1)~=length(A.data)
        wbdie('QSpace init size inconsistency');
     end

  elseif isnumeric(varargin{1})
     if isequal(varargin{end},'identity')
        [A,varargout{:}]=initQSpaceA0(varargin{1:end-1});
     elseif isequal(varargin{end},'operator')
        [A,varargout{:}]=initQSpaceH0(varargin{1:end-1});
     elseif isequal(varargin{end},'map')
        [A,varargout{:}]=initQMap(varargin{1:end-1});
     elseif isequal(varargin{end},'nosym')
        if nargin~=4, wbdie('invalid usage #8'); end
        if ~isempty(varargin{3})
             A=get_struct(...
               repmat(varargin(1),1,varargin{2}), varargin(3), []);
        else A=get_struct(); end
     elseif nargin==3 && diff(size(varargin{2}))==0
        A=initQSpaceMD(varargin{:});
     else
        A=initQSpace(varargin{:});
     end

  elseif isa(varargin{1},'QSpace') || isstruct(varargin{1})

     if isequal(varargin{end},'unity')
        A=initQSpaceUnity(varargin{1:end-1});
     else
        wbdie('invalid QSpace constructor set');
     end
  else
     wbdie('invalid QSpace constructor set');
  end

  A=struct_to_QSpace(A);
  check_QSpace(A);

end

% -------------------------------------------------------------------- %
% outsourced from usage #1 // Wb,Jan19,20

function A=struct_to_QSpace(A)

  if ~isfield(A,'Q') || ~isfield(A,'data')
     wbdie('invalid QSpace input structure !?');
  end

  if ~isfield(A,'info') && numel(A), A(1).info=[]; end
  A=class(A,'QSpace');

end

% -------------------------------------------------------------------- %
% outsourced from end of main usage above // Wb,Jan19,20

function check_QSpace(A)

  nA=numel(A);
  for k=1:nA, [q,s]=isQSpace(A(k));
     if ~q, A(k), 
        if nA>1, s=[s, sprintf(', %d/%d',k,nA)]; end
        wbdie('invalid QSpace (%s)',s);
     end
  end

end

% -------------------------------------------------------------------- %

function A = initQSpace(varargin)
% initQSpace - initialize QSpace from partial data
%
% Usage: A = initQSpace(q1, data1, q2, data2, ...)
% Example: rank-2 object with two abelian quantum numbers (QDIM=2)
%
%     A=initQSpace(...
%       [-1  0; 0  1],  1, ...
%       [ 0 -1; 1  0], -1, ...
%     );
%
% Wb,Aug08,06

  if nargin<2 || mod(nargin,2)~=0
     eval(['help ' mfilename]);
     wbdie('invalid number of arguments (usage #2)');
  end

  for i=1:2:nargin
     if ~isequal(size(varargin(i)), size(varargin(1)))
        eval(['help ' mfilename]);
        wbdie('size mismatch in Q data (usage #2)');
     end
     if ~isnumeric(varargin{i}) || ~isnumeric(varargin{i+1})
        eval(['help ' mfilename]);
        wbdie('invalid data (must be numeric, usage #2)');
     end
  end

  [rank,qdim]=size(varargin{1});
  n=nargin/2;

  qq = permute(reshape( cat(2,varargin{1:2:end}), ...
              [rank,qdim,n]), ...
       [3 2 1]);
  Q = cell(1,rank); for i=1:rank, Q{i}=qq(:,:,i); end

  A=get_struct(Q, varargin(2:2:end)', []);

end

% -------------------------------------------------------------------- %

function A = initQSpaceMD(varargin)
% initQSpaceMD - initialize QSpace from matrix
%
% Usage: A = QSpace(Q1, M, DB)
%
%    Matrix M must be square with consecutive block
%    dimensions given by DB
%
% Wb,Sep08,06

  if nargin~=3
     eval(['help ' mfilename]);
     wbdie('invalid number of arguments (usage #3)');
  end

  Q=varargin{1};
  M=varargin{2};
  D=varargin{3}; n=length(D);

  if size(Q,1)~=n || any(size(M)~=sum(D))
     eval(['help ' mfilename]);
     wbdie('dimension mismatch of arguments (usage #3)');
  end

  M = mat2cell(M, D, D);
  mark=zeros(n,n);

  for i=1:n
  for j=1:n, mark(i,j) = any(M{i,j}(:)); end, end

  if any(diag(mark)), mark=mark+eye(size(mark)); end

  [I,J]=find(mark);

  A=get_struct({Q(I,:), Q(J,:)}, M(find(mark)), []);

end

% -------------------------------------------------------------------- %

function [A,K] = initQSpaceH0(varargin)
% initQSpaceH0 - initialize QSpace from matrix
%
% Usage: [A [,K]] = QSpace(Q, M)
%
%    Matrix M must be square with quantum numbers
%    of every index specified in Q
%
% Wb,Sep08,06

  dflag=[];

  if nargin==3
     switch varargin{end}
       case 'diag', dflag=1; varargin(end)=[];
       case 'econ', dflag=0; varargin(end)=[];
       otherwise e=1;
     end
  end

  if length(varargin)~=2
     wbdie('invalid operator initialization');
  elseif size(varargin{1},1)~=size(varargin{2},1)
     eval(['help ' mfilename]);
     wbdie('dimension mismatch of arguments (usage #4)');
  end

  [Q,K,D]=uniquerows(varargin{1}); n=length(K);
  mark=zeros(n,n); M=cell(n,n);

  for i=1:n, for j=1:n
     a=varargin{2}(K{i},K{j});
     if any(a(:))
        M{i,j}=full(a); mark(i,j)=1;
     end
  end, end

  for g=1:length(K)
     p=sort(reshape(K{g},1,prod(size(K{g})))); if ~isequal(p,K{g})
       wblog('WRN','index sets not sorted in uniquerows() ???');
       K{g}=p;
     end
  end

  if isempty(dflag) && any(diag(mark)) || ~isempty(dflag) && dflag
     for i=1:n
        if mark(i,i), continue; end
        M{i,i}=full(varargin{2}(K{i},K{i}));
        mark(i,i)=1;
     end
  end

  [I,J]=find(mark);

  A=get_struct({Q(I,:), Q(J,:)}, M(find(mark)), []);

end

% -------------------------------------------------------------------- %

function varargout = initQMap(varargin)
% initQMap - initialize QSpace from matrix
%
% Usage: [A [,K]] = QSpace(Q1,Q2,...,M)
%
%    M is regular matrix
%    with quantum numbers of every index (i,j) specified in Q1,Q2
%
% generalized to arbitrary rank>=2 // Wb,Mar25,14
% Wb,Aug06,07

  n=nargin; if n<3, e=-100;
  else e=0; M=varargin{end};
     for i=1:n
        if ~isnumeric(varargin{i}), e=i;
        elseif i<n && size(varargin{i},1)~=size(M,i), e=100*i; end
     end
  end
  if e
     eval(['help ' mfilename]);
     wbdie('dimension mismatch of arguments (initQMap [e=%g])',e);
  end

  n=n-1; s=zeros(1,n);
  for i=1:n
    [Q{i},K{i},D{i}]=uniquerows(varargin{i});
    s(i)=size(Q{i},1);
    for j=1:s(i)
       p=sort(reshape(K{i}{j},1,[])); if ~isequal(p,K{i}{j})
         wblog('WRN','index sets not sorted in uniquerows() ???');
         K{i}{j}=p;
       end
    end
  end

  mark=zeros(s); data=cell(s); i=cell(1,n); I=cell(1,n);
  if n>2, M=full(M); end

  for l=1:prod(s), [i{:}]=ind2sub(s,l);
     for j=1:n, I{j}=K{j}{i{j}}; end
     a=M(I{:});
     if any(a(:))
        data{i{:}}=full(a); mark(i{:})=1;
     end
  end

  l=find(mark); [i{:}]=ind2sub(s,l);
  for j=1:n, I{j}=Q{j}(i{j},:); end

  A=get_struct(I,data(l),[]);

  if nargout<2
       varargout={A};
  else varargout={A,K{:}}; end

end

% -------------------------------------------------------------------- %
% build identity operator from other operator basis

function E = initQSpaceUnity(A)

  if ~nargin || ~isa(A,'QSpace') && (~isstruct(A) || ~isfield(A,'Q'))
     wbdie('invalid usage'); end
  q=numel(A); if q~=1
     wbdie('invalid usage (single rank-2 object expected; n=%d)',q); end
  q=numel(A.Q); if q~=2
     wbdie('invalid usage (rank-2 object expected; r=%d)',q); end

  [q1,d1]=getQDimQS(A,1);
  [q2,d2]=getQDimQS(A,2);

  QD=[ q1, d1; q2, d2 ];
  [qd,I,d]=uniquerows(QD(:,1:end-1)); n=length(I);

  Q=qd; D=zeros(n,1);

  for i=1:n, di=QD(I{i},end); 
     if d(i)>1 && any(diff(di)), di
        wbdie('invalid rank-2 object / severe QSpace inconsistency'); end
     D(i)=di(1);
  end

  E=get_struct();
  E.Q={Q,Q};
  E.data=cell(1,n); for i=1:n, E.data{i}=eye(D(i)); end

end

% -------------------------------------------------------------------- %

function [A,Aloc,I1,I2] = initQSpaceA0(varargin)

% init MPS block A0 representing identity (returns LRs order)
% Usage: A0= initQSpace(QL,Qs1,Qs2,...);
%
%   The tensor product of QL and possibly multiple Qs spans
%   the full Hilbert space; they may not be unique yet; grouping
%   them is used to determine their degeneracies.
%
% NB! returns (QL,QR=Qtot,Qs1,Qs2,...) order
%
% Wb,Sep10,06

  if nargin<2, wbdie('invalid usage'); end

  getopt('init',varargin);
     Rlast=getopt('-Rlast');
     mergeloc=getopt('mergeloc');
  varargin=getopt('get_remaining');

  nqin=length(varargin);
  for i=1:nqin
     if ~isnumeric(varargin{i})
     wbdie('remaining non Q-matrix type of argument !?'); end
  end

  QDIM=size(varargin{1},2);
  for i=2:nqin
     if size(varargin{i},2)~=QDIM
        wbdie('QDIM mismatch (%d,%d)',QDIM,size(varargin{i},2));
     end
  end

  if mergeloc
     for i=2:nqin
     ee{i-1}=ones(size(varargin{i},1),1); end

     for i=2:nqin
     ei=ee; ei{i-1}=varargin{i}; QQ{i-1}=mkron(ei{:}); end

     varargin={ varargin{1}, cellarrsum(QQ) };
     nqin=length(varargin);
  end

  dd=zeros(1,nqin); cn=cell(1,nqin); ee=cn; eu=cn; du=cn;
  for i=1:nqin
      dd(i)=size(varargin{i},1);
      ee{i}=ones(dd(i),1);
      [~,~,du{i}]=uniquerows(varargin{i});
      du{i}=double(reshape(du{i},[],1));
      eu{i}=ones(length(du{i}),1);
  end

  QQ=cn; DD=cn;
  for i=1:nqin
      ei=ee; ei{i}=varargin{i}; QQ{i}=mkron(ei{:});
      ei=eu; ei{i}=du{i};       DD{i}=mkron(ei{:},'rowmajor');
  end

  DD=cat(2,DD{:});

  [D,QDIM]=size(QQ{1}); E=speye(D,D);

  [Q1,I1,D1]=uniquerows(cat(2,QQ{:}));   m=length(D1);
  [Q2,I2,D2]=uniquerows(cellarrsum(QQ)); n=length(D2);

  M=cell(m,n); mark=zeros(m,n);

  for i=1:m
     if D1(i)~=prod(DD(i,:))
        wbdie('failed to determine dimensions of Q [%d; %s]', ...
        D1(i), vec2str(DD(i,:)));
     end

     for j=1:n
        a=E(I1{i},I2{j}); if any(a(:))
           mark(i,j)=1;
           M{i,j}=reshape(full(a), [ DD(i,:), D2(j)]);
        end
     end
  end

  [I,J]=find(mark); IJ=sub2ind(size(mark),I,J);

  Q1=Q1(I,:);
  Q1=mat2cell(Q1, size(Q1,1), QDIM(ones(1,nqin)));

  A=get_struct({Q1{:}, Q2(J,:)}, M(IJ), []);

  if ~Rlast
     p=1:(nqin+1); p=p([1 end 2:end-1]);
     A=permuteQS(A,p);
  end

  if nargout>1 
     if length(varargin)==2
          Aloc=QSpace(varargin{2},eye(size(varargin{2},1)),'operator');
     else Aloc=QSpace(varargin{2:end},'-Rlast','identity'); end
  end

end

% -------------------------------------------------------------------- %

