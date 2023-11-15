function [i,bfac,e]=sameup2fac(A,B,varargin)
% function [i,bfac,e]=sameup2fac(A,B [,eps])
%
%   checks whether A = bfac * B
%   within numerical noise (default: eps=1e-12).
%
% Options
%
%   '-ev'  check whether `eigenvectors' A and B span the same state space,
%          i.e., A = B * bfac with bfac a diagonal matrix here.
% 
% Wb,Jun29,10

  if nargin<2
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  getopt('init',varargin);
     evflag=getopt('-ev');
  eps=getopt('get_last',1E-10);

  if ~isnumeric(A) || ~isnumeric(B)
     wbdie('invalid input (numeric data expected)'); end

  i=0; e=nan;
  if ~isequal(size(A),size(B)), return; end

  if evflag
     B2=B'*B; bfac=inv(B2)*(B'*A);

     e=norm(bfac-diag(diag(bfac)))/norm(bfac);
     if e>eps, return; end; e=norm(A-B*bfac)/norm(B2);
     if e>eps, return; end

  else
     a=A(:); b=B(:);
     b2=b'*b;

     bfac=full((b'*a)/b2);

     if b2==0, wbdie('got |B|=0');
     else
        a2=a'*a; q=b2/a2;
        if q<1E-12, wblog('WRN','got b2/a2 = %.3g',q); end
     end

     if ~isreal(bfac), q=[real(bfac), imag(bfac)];
        if     abs(q(2)/q(1))<1E-14, bfac=q(1);
        end
     end

     e=norm(a-bfac*b)/sqrt(b2);
     if e<1E-8 && e>1E-12, wblog('WRN','got small deviation (%.3g)',e); end
     if e>eps, return; end
  end

  i=1;

end

