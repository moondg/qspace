function x=trace(A,I,varargin)
% Function x=trace(A,I [,opts])
%
%   trace of QSpace operator
%
% Usage 1: A=trace(A)   - regular trace of rank-2 object
% Usage 2: A=trace(A,I)
%
%    generalized trace with pairwise specification of indizes
%    to contract I = [i1 i2; j1 j2; ... ]
%
% Wb,Sep11,06 ; Wb,Apr10,15

   if ~nargin
      eval(['help ' mfilename]);
      if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
   end

   getopt('init',varargin);
      Qflag=getopt('-Q');
   getopt('check_error');

   if isscalar(A)
      if nargin>1 && ~isempty(I)
      error('Wb:ERR','cannot contract scalar'); end
      x=A.data{1}; return
   elseif isempty(A.Q)
      if Qflag, x=QSpace({},{0}); else x=0; end
      return
   end

   r=rank(A); cgflag=gotCGS(A); isd=isdiag(A);

   if nargin<2, x=0;
      if mod(r,2), error('Wb:ERR',...
      '%s requires even-rank object (%x).',mfilename,r); end

      Q1=cat(2,A.Q{1:r/2});
      Q2=cat(2,A.Q{r/2+1:r});

      for i=1:length(A.data)
         if isequal(Q1(i,:), Q2(i,:)), d=A.data{i}; s=size(d); q=numel(s);
            if q>2, if mod(q,2), s(end+1)=1; q=q+1; end
                q=q/2; if ~isequal(s(1:q),s(q+1:end))
                error('Wb:ERR','operator dimensions must be symmetric'); end
                d=reshape(d,prod(s(1:q)),[]);
            end
            if isd>1
                 d=sum(d);
            else d=trace(d); end

            if cgflag
               cg=A.info.cgr(i,:); m=numel(cg);
               for j=1:m
                  if ~isempty(cg(j).cgw)
                     q=mpfr2dec(cg(j).cgw) .* mpfr2dec(cg(j).cgt);
                     d=d*sum(q(:));
                  elseif ~isempty(cg(j).type) || ~isempty(cg(j).qset)
                     error('Wb:ERR','\n   ERR invalid CGR_ABELIAN');
                  end
               end
            end
            x=x+d;
         end
      end

      if Qflag, x=QSpace({},{x}); end
      return
   end

 % generalized trace
 % pairwise specification of indizes to contract I = [i1 i2; j1 j2; ... ]
 % => no longer applicable (need mex file)
 % see Archive/trace_141112.m for old m-script implementation
 % Wb,Nov12,14

   error('Wb:ERR','\n   ERR invalid usage');
end

