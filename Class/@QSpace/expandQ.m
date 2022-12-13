function qq=expandQ(A,k)
% Function qq=expandQ(A,k)
% deprecated: see getQRange.m

% [Q,I,d]=uniquerows(A.Q{k});
% qq=I; n=length(I);
% for i=1:n, qq{i}=Q(ones(size(A.data{I{i}(1)},1),1),:); end

  [Q,d]=getQDimQS(A,k); % k=='op' is allowed (consider both dimenions)
  n=length(d); qq=cell(n,1);

  for i=1:n, qq{i}=Q(i(ones(d(i),1)),:); end

  qq=cat(1,qq{:});

end

