function A=symperm(A,perm)
% function A=symperm(A,perm)
%
%    Permute symmetries according to specified permutation.
%
% Wb,Feb06,14

  if nargin~=2
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  ns=numel(perm);
  if size(perm,1)~=1 || ~isequal(sort(perm(:)'),1:ns), perm
     error('Wb:ERR','\n   ERR got invalid permutation');
  end

  for k=1:numel(A)
     [d,s]=getsym(A(k),'-d');
     if numel(s)~=ns, error('Wb:ERR',...
       '\n   ERR got invalid permutation (%g/%g)',numel(s),ns);
     end

     if ~isempty(A(k).info.qtype)
        s=s(perm)'; s(2,1:end-1)={','};
        A(k).info.qtype=[s{:}];
     end

     if ~isempty(A(k).info.cgr)
        A(k).info.cgr=A(k).info.cgr(:,perm);
     end

     if ~isempty(A(k).Q) && ~isempty(A(k).Q{1}) 
        n=size(A(k).Q{1},2); if n~=sum(d), error('Wb:ERR',...
          '\n   ERR inconsistent Q data (%g/%g)',n,size(A(k).Q{1},2));
        end
        p=mat2cell(1:n,1,d); p=[p{perm}];
        for j=1:numel(A(k).Q)
           A(k).Q{j}=A(k).Q{j}(:,p);
        end
     end
  end

end

