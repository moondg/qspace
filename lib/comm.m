function C=comm(A,B)
% Function: C=comm(A,B)
%
%    calculate commuator [A,B]=AB-BA
%    of two objects A and B
%
% Wb,Jun19,07

   if nargin<2, B=A'; end
   C=A*B-B*A;

end

