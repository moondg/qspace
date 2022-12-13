function x=norm(A,varargin)
% function x=norm(A)
% Wb,Dec05,11

   if ~isa(A,'SymOp'), error('Wb:ERR','\n   ERR invalid usage'); end

   getopt('init',varargin);
      odiag=getopt('-offdiag');
   getopt('check_error');

   x=zeros(size(A)); n=numel(A);
   for i=1:n
      q=A(i).op; if odiag, q=q-diag(diag(q)); end
      x(i)=norm(q,'fro');
   end

end

