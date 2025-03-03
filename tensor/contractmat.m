function C=contractmat (A,B,ia,P_);
% function C = contractmat (A,B,ia [,pout])
%
%     Contract single index ia of tensor A with matrix B
%     while maintaining the index order in A. The resulting C
%     can be optionally permuted by the trailing input argument pout.
%
% See contract() for more general contraction of two tensors.
% F.Verstraete / adapted by AWb

  sa=size(A); ra=max(numel(sa),ia); sa(end+1:ra)=1;
  sb=size(B); q=[sa(ia), sb(1)];
     if numel(sb)~=2, wbdie('invalid usage (got non-matrix B)'); end
     if diff(q), wbdie('dimension mismatch (%d/%d)',q); end
  Ia=1:ra; Ia(ia)=[]; P=[Ia ia]; iP(P)=1:ra;

  AA=permute(A,P);
  C=reshape(AA,[prod(sa(Ia)) sa(ia)])*B;
  C=reshape(C, [sa(Ia)       sb(2 )]);

  if nargin==4, iP=iP(P_); end

  C=permute(C,iP);

end

