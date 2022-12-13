function [a,I,D]=groupvec(A)
% function [a,I,D]=groupvec(A)
%
%    replica of mex-function uniquerows()
%    for non-numeric data type such as graphics handles.
%
% Wb,Aug03,16

  n=numel(A); if ~n, a=A; I={}; D=[]; return; end

  [aa,is]=sort(A(:));

  mark=ones(n,1);
  for i=2:n
     if aa(i-1)==aa(i), mark(i)=0; end
  end

  ix=find(mark); a=aa(ix);
  m=numel(ix); I=cell(1,m); D=zeros(1,m);

  for i=1:m
     if i<m
          I{i}=is(ix(i):ix(i+1)-1);
     else I{i}=is(ix(i):n);
     end
     D(i)=numel(I{i});
  end

end

