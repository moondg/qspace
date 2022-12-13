function Q=contractmat (A, B, ma, permQ);
% Function - contractmat
% Usage: Q = contractmat ( A, B, ma, permQ );
%
%     contract single index ma of tensor A with matrix B and
%     replace this index ma with the one labeling the columns of B
%
% AWb, F.Verstraete

% Wb - see contract() for the more general
% version of contracting two arbitrary tensors

  sa=size(A);na=max(size(sa));Ia=[1:na];Ia(ma)=[];
  sb=size(B);
  AA=permute(A,[Ia ma]);
  Q=reshape(AA,[prod(sa(Ia)) sa(ma)])*B;
  Q=reshape(Q,[sa(Ia) sb(2)]);

  Q=permute(Q,[1:ma-1, na, ma:na-1]);

  if nargin==4, Q=permute(Q,permQ); end

end

