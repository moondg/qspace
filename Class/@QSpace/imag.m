function A=imag(A)
% Function: A=imag(A)
%
%    extract imaginary data of operator in QSpace format
%
% Wb,Oct19,11

  if nargin<1 || nargout>1
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

  for k=1:numel(A)
     n=length(A(k).data);
     for i=1:n, A(k).data{i}=imag(A(k).data{i}); end
  end

end

