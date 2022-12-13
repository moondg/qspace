function [s,dd]=isscalarop(A,varargin)
% function [s,dd]=isscalarop(A [,opts])
%
%    requirements for scalar operators:
%     (1) rank-2
%     (2) block-diagonal in symmetries
%     (3) all cgc spaces are identities
%
% Options
%
%   -l   lenient (also accept rank-3 operators)
%
% See also @QSpace/issingleop for non-scalar, yet still rank-2 operators.
% Wb,Jul05,12

  getopt('init',varargin);
     lflag=getopt('-l');
  getopt('check_error');

  if numel(A.data)==0 && numel(A.Q)==0, s=1; return ; end

  r=numel(A.Q);
  if lflag
     if r<2 || r>3 || norm(A.Q{1}-A.Q{2})>1E-12
        s=0; return
     end
     if r==3 && norm(diff(A.Q{3},[],1))>1E-12
        s=0; return
     end
  else
     if r~=2 || norm(A.Q{1}-A.Q{2})>1E-12 || dim(A,3,'-f')~=1
        s=0; return
     end
  end

  if ~gotCGS(A), s=1; return; end

  n=numel(A.info.cgr);

  for i=1:n, c=A.info.cgr(i);
     if mod(numel(c.size),2), s=0; return; end
     if norm(diff(reshape(double(c.size),[],2),[],2)), s=0; return; end
     if norm(diff(reshape(c.qdir,[],2),[],2)-2), s=0; return; end
     if norm(diff(reshape(c.qset,[],2),[],2)), s=0; return; end
  end

  s=1; if nargout<2, return; end

  n=numel(A.data); dd=nan(size(A.data));

  for i=1:n, q=A.data{i};
      e=norm(q-q(1)*eye(size(q)));
      if e<1E-12, dd(i)=q(1); end
  end

  if none(isnan(dd)), s=2; end

end

