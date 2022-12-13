function A=conj(A)
% Function: A=conj(A)
%
%    get complex conjugate of operator in QSpace format;
%    overloading the conj routine
%
% Wb,Jun02,08

  if nargin<1 || nargout>1
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

  for k=1:numel(A), A(k)=SetCC(A(k)); end

end

