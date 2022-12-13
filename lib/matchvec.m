function [ia,ib,I]=matchvec(A,B,sflag)
% function [ia,ib,I]=matchvec(A,B [,'-s'])
%
%    replica of mex-function matchIndex()
%    for non-numeric data type such as graphics handles.
% 
% The returned structure I contains the fields
%
%    ix1 
%    ix2   index for the unmatched  entries in A and B, respectively.
%    m     whether multiplicity (non-uniqueness of entries was encountered).
%
% Wb,Aug03,16

  if nargin>2
     if nargin>3 || ~isequal(sflag,'-s')
        wberr('invalid usage'); end
     sflag=1;
  else sflag=0;
  end

  n=[numel(A), numel(B)];
  if any(n==0)
     ia=[]; ib=[]; I.ix1=1:n(1); I.ix2=1:n(2); I.m=0;
     return
  end

  [aa,Ia,Da]=groupvec(A); na=numel(aa);
  [bb,Ib,Db]=groupvec(B); nb=numel(bb);

  i=1; j=1; l=1; nab=max(na,nb); m=zeros(nab,1);
  IJ=cell(nab,1); xa=zeros(na,1); xb=zeros(nb,1);

  while i<=na && j<=nb
     if     aa(i)<bb(j), i=i+1;
     elseif aa(i)>bb(j), j=j+1;
     else
        IJ{l}=vkron(Ia{i},Ib{j}); m(l)=size(IJ{l},1);
        xa(i)=j; xb(j)=i;  i=i+1; j=j+1; l=l+1;
     end
  end

  if nargout>2
     i=find(xa==0); I.ix1=sort(cat(1,Ia{i}))';
     i=find(xb==0); I.ix2=sort(cat(1,Ib{i}))';
     I.m=sum(m(1:l-1)-1);
  end

  q=cat(1,IJ{:}); if sflag, q=sortrows(q); end
  if ~isempty(q)
       ia=q(:,1)'; ib=q(:,2)';
  else ia=[]; ib=[]; end

end

