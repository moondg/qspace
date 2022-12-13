function A=chopd(A,n)
% function A=chopd(A [,n])
%
%    Chop noise from double precision data
%    by converting to single and back to double.
%    NB! used e.g. to sort and group float data.
%
%    If second argument n is specified, A is rounded to n digits.
%
% See also chop.m
% Wb,Nov05,08

% NB! for complex numbers, rather use chop() than *this.

 % skip small imaginary component if present
   if ~isreal(A)
      a=[norm(real(A(:))), norm(imag(A(:)))];
      if a(2)<1E-9*a(1), A=real(A); end
   end

   if nargin<2, one=1; if ~isreal(A), one=complex(one,one); end
      i=find(A>0); if ~isempty(i), A(i)=double(single(A(i)+one))-one; end
      i=find(A<0); if ~isempty(i), A(i)=double(single(A(i)-one))+one; end
   elseif ischar(n)
      if ~isequal(n,'-i'), wbdie('invalid usage (option %s)',n); end
      a=max(abs(A(:))); if a<1, a=1; end
      A_=round(A);
      i=find(abs(A-A_) < 1E-8*a); if ~isempty(i), A(i)=A_(i); end
   else
      fac=10^n;
      A=round(A*fac)/fac;
   end

end

