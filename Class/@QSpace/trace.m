function x=trace(A,I,varargin)
% function x=trace(A,I [,opts])
%
%   trace of QSpace tensor
%
% Usage 1: A=trace(A)   - regular trace of rank-2 object
% Usage 2: A=trace(A,I [,J])
%
%    generalized trace with pairwise specification of indizes
%    to contract I = [i1 i2; j1 j2; ... ]
%    or, if J is present, I with J.
%
% Wb,Sep11,06 ; Wb,Apr10,15 ; Wb,Jul16,25

   if ~nargin
      if ~helpthis(nargout,varargin{:}), wbdie('invalid usage'); end
     return
   end

   getopt('init',varargin);
      Qflag=getopt('-Q');
   J=getopt('get_last',[]);

   n=numel(A);
   if n~=1, wbdie('invalid usage (%d entries in A)',n); end

   if isscalar(A)
      if nargin>1 && (~isempty(I) || ~isempty(J))
         wbdie('cannot contract scalar'); end
      if Qflag, x=A; else x=A.data{1}; end
      return
   elseif isempty(A.Q)
      if Qflag, x=QSpace({},{0}); else x=0; end
      return
   end

   r=rank(A); cgflag=gotCGS(A); isd=(isdiag(A)>1);

   if nargin<2, x=0;
      if mod(r,2)
         wbdie('%s requires even-rank object (%d)',mfilename,r); end

      if r>2
         d=getDimQS(A); d=d(end,:);
         if any(d(1:2)~=1) && all(d(3:end)==1), i1=1; i2=2;
         else
            i1=1:2:r; i2=2:2:r;
            if ~isequal(d(i1),d(i2)), i1=1:r/2; i2=r/2+1:r; end
         end
      else i1=1; i2=2;
      end

      Q1=cat(2,A.Q{i1});
      Q2=cat(2,A.Q{i2});

      for i=1:length(A.data)
         if isequal(Q1(i,:), Q2(i,:)), d=A.data{i}; s=size(d); q=numel(s);
            if q>2, if mod(q,2), s(end+1)=1; q=q+1; end
                q=q/2; if ~isequal(s(1:q),s(q+1:end))
                wbdie('operator dimensions must be symmetric'); end
                d=reshape(d,prod(s(1:q)),[]);
            end
            if isd
                 d=sum(d);
            else d=trace(d); end

            if cgflag
               cg=A.info.cgr(i,:); m=numel(cg);
               for j=1:m
                  if ~isempty(cg(j).cgw)
                     q=cg(j).cgw .* cg(j).cgt;
                     d=d*sum(q(:));
                  elseif ~isempty(cg(j).type) || ~isempty(cg(j).qset)
                     wbdie('invalid CGR_ABELIAN');
                  end
               end
            end
            x=x+d;
         end
      end

      if Qflag, x=QSpace({},{x}); end
      return
   end

   if nargin==2
      if size(I,2)~=-2, wbdie('invalid usage'); end
      J=I(:,2); I=I(:,1);
   elseif ~isequal(size(I),size(J))
      wbdie('invalid usage (I and J of different size)');
   end

   x=traceQS(A,I,J);
   x=QSpace(x);
end

