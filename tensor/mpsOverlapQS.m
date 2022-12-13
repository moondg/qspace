function x = mpsOverlapQS(A,B);
% function x = mpsOverlapQS(A,B);
%
%    calculate MPS overlap <A|B> in between the two MPS A and B.
%    Convention: LRs index order in A-tensors.
%
% Wb,Aug05,06 ; Wb,Apr08,14

  if nargin==1, B=A; end

  if ~nargin || ~isQSpace(A) || nargin>1 && ~isQSpace(B)
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end
  if ~isvector(A) || ~isvector(B) || numel(A)~=numel(B)
     error('Wb:ERR','\n   ERR incompatible MPS states');
  end

  N=numel(A);

  d1={ getDimQS(A(1)), getDimQS(B(1)) };
  d2={ getDimQS(A(N)), getDimQS(B(N)) };

  if size(d1{1},2)~=3 || size(d1{2},2)~=3
     error('Wb:ERR',['\n   ERR invalid MPS state(s) (got rank ' ... 
    '%g/%g tensors (k=1) !??)'], size(d1{1},2), size(d1{2},2));
  end
  if size(d2{1},2)~=3 || size(d2{2},2)~=3
     error('Wb:ERR',['\n   ERR invalid MPS state(s) (got rank ' ... 
    '%g/%g tensors (k=N) !??)'], size(d2{1},2), size(d2{2},2));
  end

  if d1{1}(end,1)~=1 || d1{2}(end,1)~=1
     s=sprintf('got non-scalar L-vac state (%g,%g) !??',...
        d1{1}(end,1),d1{2}(end,1));
     if d1{1}(1,1)~=1 || d1{2}(1,1)~=1 || ~isequal(d1{1}(:,1),d1{2}(:,1))
        error('Wb:ERR','\n   ERR %s',s); 
     else wblog('WRN','%s',s); end
  end

  if d2{1}(end,2)~=1 || d2{2}(end,2)~=1
     s=sprintf('got non-scalar R-vac state (%g,%g) !??',...
        d2{1}(end,2),d2{2}(end,2));
     if d2{1}(1,2)~=1 || d2{2}(1,2)~=1 || ~isequal(d2{1}(:,2),d2{2}(:,2))
          error('Wb:ERR','\n   ERR %s',s); 
     else wblog('WRN','%s',s); end
  end

  if ~isvector(A) || ~isvector(B) || numel(A)~=numel(B)
     error('Wb:ERR','\n   ERR incompatible MPS states');
  end

  k=1; x=contractQS(A(k),'13*',B(k),'13');

  for k=2:N-1
     x=contractQS(A(k),'13*',contractQS(x,2,B(k),1),'13');
     if isempty(x.data), x=0; return; end
  end

  k=N; x=contractQS(A(k),'123*', contractQS(x,2,B(k),1),'123');

  x=getscalar(QSpace(x));

end

