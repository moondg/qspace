function C=acomm(A,B)
% Function: C=acomm(A,B)
%
%    calculate anti-commuator {A,B}=AB+BA
%    of two objects A and B
%
% Wb,Jun19,07

   if nargin<2, B=A'; end
   C=A*B+B*A;

end

