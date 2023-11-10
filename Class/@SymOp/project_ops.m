function [A,x]=project_ops(A,B,eps)
% function [A,x]=project_ops(A,B [,eps])
%     A = A - x*|B> with x = <B|A>/|B|^2
% Wb,Sep27,12

% NB! renamed function from project() to project_ops()
% since project() collided with license error during mcc compilation:
%    No licenses available for toolbox file // !?!?#_&!_*#$&)*(_!?!?
%   '/software/opt/xenial/x86_64/matlab/2016a/toolbox/map/mapdisp/project.m'.
% Compilation terminated.
% Wb,Sep27,17

  if nargin<2 || isa(A,'SymOp') && numel(A)~=1 || isa(B,'SymOp') && numel(B)~=1
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end
  if nargin<3, eps=0; end

  if isa(A,'SymOp')
     if isa(B,'SymOp'), B=B.op; end
     x=olap(A.op,B,'-nB');
     if abs(x)>eps, A.op=A.op-x*B; A.istr=''; else x=0; end
  else
     if isa(A,'SymOp'), A=A.op; end
     x=olap(A,B.op,'-nB');
     if abs(x)>eps, B.op=A-x*B.op; B.istr=''; else x=0; end
     A=B;
  end

end

