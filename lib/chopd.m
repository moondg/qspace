function A=chopd(A,n)
% function A=chopd(A [,n])
%
%    Chop noise from double precision data in A (which may also be a
%    cell array of such data) by converting to single and back to double.
%    May be used, e.g., to sort and group floating-point data.
%    If second argument n is specified, A is rounded to n digits.
%
% See also chop.m
% Wb,Nov05,08

% NB! for complex numbers, rather use chop() than *this.

  if nargin<2, n=[];
  elseif ischar(n)
         if ~isequal(n,'-i'), n='i'; % fix integers ('-i') // Wb,Oct21,22
     elseif ~isequal(n,'-c'), n='c';
     else wbdie('invalid usage (option %s)',n); end
  elseif ~isnumeric(n) || numel(n)~=1
     wbdie('invalid usage'); 
  end

  if ~iscell(A)
      A=chopd_1(A,n);
  else
     for i=1:numel(A), A{i}=chopd_1(A{i},n); end
     if ~nargout, A{:}, clear A, end
  end

end

% -------------------------------------------------------------------- %
% Wb,Nov05,08

function A=chopd_1(A,n)

   if ~isreal(A)
      a=[norm(real(A(:))), norm(imag(A(:)))];
      if a(2)<1E-9*a(1), A=real(A); end
   end

   if isempty(n), a=1; if ~isreal(A), a=complex(a,a); end
      i=find(A>0); if ~isempty(i), A(i)=double(single(A(i)+a))-a; end
      i=find(A<0); if ~isempty(i), A(i)=double(single(A(i)-a))+a; end
   elseif isequal(n,'i')
      a=max(abs(A(:))); if a<1, a=1; end
      A_=round(A);
      i=find(abs(A-A_) < 1E-8*a); if ~isempty(i), A(i)=A_(i); end
   elseif ~isequal(n,'c')
      fac=10^n;
      A=round(A*fac)/fac;
   end

end

% -------------------------------------------------------------------- %

