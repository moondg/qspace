function [i,I] = hasQOverlap(A,i1,B,i2,flag)
% function: [ss,qq] = hasQOverlap(A,i1,B,i2)
% check overlap in QIDX
% Wb,Jun22,07

  Q1=uniquerows(A.Q{i1});
  Q2=uniquerows(B.Q{i2}); i=0;

  [ia,ib,I]=matchIndex(Q1,Q2);

  I.ia=ia;
  I.ib=ib;

  n=length(fieldnames(I));
  I=orderfields(I,[n-1,n,1:n-2]);

  if exist('flag','var')
     if     isequal(flag,'>'), if isempty(I.ix2) i=1; end
     elseif isequal(flag,'<'), if isempty(I.ix1) i=1; end
     else   if ~isempty(ia), i=1; end
     end
  end

end

