function x=MPS_overlap(A,B)
% function x=MPS_overlap([A,] B)
%
%    Compute MPS overlap <A|B>.
%    If only a single argument (MPS) is specified
%    by default A=B, i.e., this computes norm2(A) = <A|A>.
%
% Wb,Feb11,25

% adpated from MPO_overlap // Wb,Feb11,25

  s=size(A); L=prod(s); ic='!1*';
  if s(1)~=1 || numel(s)>2
     wbdie('unexpected MPS QSpace array dimensions'); end

  if nargin<2, B=A; % compute norm2(A) = <A|A> // tags: MPS_norm2()
  elseif ~isequal(s,size(B)), wbdie('MPS size mismatch'); end

  for k=L:-1:1,  if k==1, ic='*'; end
      Q=B(k); if k<L, Q={Q,X}; end
      X=contract(A(k),ic,Q);
  end
  x=getscalar(X);

end

