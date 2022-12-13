function [ok,ii]=checkdim(A,pp)

  ok=1; ii=[];

% get size
  r=length(A.Q); m=length(A.data); sa=zeros(m,r);
  for i=1:m
  s=size(A.data{i}); sa(i,1:length(s))=s; end

  if nargin==1, pp=1:r; end
  for p=pp
     [Q,I,D]=uniquerows(A.Q{p});
     for i=1:length(D)
         s=sa(I{i},p); if any(s~=s(1)), ok=0; ii=I{i}; return, end
     end
  end

return

