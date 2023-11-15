function [A,isd]=sort(A,perm)
% function [A,isd]=sort(A,perm)
%
%    sort records in QSpace according to q-labels after using
%    permutation perm (default: Id). Alternatively, perm may
%    also be specified as 
%
%      '-E' (sort by energy) or
%      '-R' (same as '-E', but sorts in descending order).
%      '-S' sort wrt. to size of data+cgr
%
% Wb,Sep12,06

  r=length(A(1).Q); if nargin<2, perm=1:r; end
  nA=numel(A); isd=0;

  for k=1:numel(A)

     if isnumeric(perm)

        [x,is]=sortrows(cat(2,A(k).Q{perm}));

     elseif isequal(perm,'-E') || isequal(perm,'-R')

        isd(k)=isdiag(A(k),'-d');
        if isd(k)~=1 && isd(k)~=2
           wblog('WRN','got non-diagonal operator => calling eig()');
           [~,I]=eigQS(A(k)); A(k)=I.EK;
           isd(k)=2;
        end

        n=numel(A(k).data); ee=zeros(n,1);
        if isd(k)==1
             for i=1:n, ee(i)=min(diag(A(k).data{i})); end
        else for i=1:n, ee(i)=min(A(k).data{i}); end
        end

        if isequal(perm,'-R'), o={'descend'}; else o={}; end
        [ee,is]=sort(ee,o{:});

     elseif isequal(perm,'-s') || isequal(perm,'-S')

        n=numel(A(k).data); ss=zeros(n,2);
        for i=1:n, ss(i,1)=sizeof(A(k).data{i}); end

        if isfield(A(k).info,'cgr')
           cgr=A(k).info.cgr; m=size(cgr,2); ss(:,2:1+m)=1;
           for i=1:n
           for j=1:m, s=double(cgr(i,j).size); l=length(s);
               if l, ss(i,2:1+l)=ss(i,2:1+l).*s; end
           end
           end
        end

        if isequal(perm,'-S'), o={'descend'}; else o={}; end
        [ss,is]=sort(prod(ss,2),o{:});

     else perm, wbdie('invalid perm');
     end

     for i=1:r, A(k).Q{i}=A(k).Q{i}(is,:); end

     A(k).data=A(k).data(is); if gotCGS(A(k))
     A(k).info.cgr=A(k).info.cgr(is,:); end

  end

end

