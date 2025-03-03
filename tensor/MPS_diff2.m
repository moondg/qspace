function dx2=MPS_diff2(A,B)
% function dx2=MPS_diff2([A,] B)
%
%    Compute MPS difference |A-B|^2 = trace((A-B)^\dagger(A-B))
%
% Wb,Feb11,25

% adpated from MPS_overlap // Wb,Feb11,25

  if nargin~=2, wbdie('invalid usage'); end
  if isequal(A,B), dx2=0; return; end

  s=size(A); L=prod(s); ic='!1*';
  if s(1)~=1 || numel(s)>2
     wbdie('unexpected MPS QSpace array dimensions'); end
  if ~isequal(s,size(B)), wbdie('MPS size mismatch'); end

  XR=QSpace(3,1);
  for k=L:-1:1,  if k==1, ic='*'; end
      Q=A(k); if k<L, Q={Q,XR(1)}; end; XR(1)=contract(A(k),ic,Q); % <A|A>
      Q=B(k); if k<L, Q={Q,XR(2)}; end; XR(2)=contract(A(k),ic,Q); % <A|B>
              if k<L, Q{2}=XR(3) ; end; XR(3)=contract(B(k),ic,Q); % <B|B>
      if rank(XR(2))>2, wbdie('got itag mismatch across A and B'); end
  end

  a2=getscalar(XR(1));
  ab=getscalar(XR(2));
  b2=getscalar(XR(3));

  dx2 = a2+b2-2*real(ab);

end

