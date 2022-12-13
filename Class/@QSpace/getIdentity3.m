function E=getIdentity3(A,idx,varargin)
% function E=getIdentity3(A,idx [,tag,opts])
%
%    Map pair of open indices to be combined as in delta_{ij}
%    into single index in E_(i,j,k) by mapping (i=j -> k).
%    The added third index therefore always transforms like a scalar
%    for all symmetry sectors present, i.e., has q-labels all zero.
%
% Options
%
%    --rho normalize E
%          such that it resembles density matrix in (i,j) for fixed k.
%    tag   tag assigned to the new trailing index.
%
% Wb,Aug28,20

  getopt('init',varargin);
     Rflag=getopt('--rho');
  tag=getopt('get_last','');

  E=QSpace(getIdentityQS(A,idx));
  E=struct(appendSingletons(E,-1)); nd=numel(E.data); 

  ed=ones(nd,2);
  for i=1:nd, ed(i,:)=[E.data{i}(1), size(E.data{i},1)]; end

  dc=[0; cumsum(ed(:,2))]; D=dc(end);
  if Rflag
     [q_,d_,dq]=getQDimQS(E,1);
     ed(:,1)=ed(:,1)./prod(dq,2);
  end

  for i=1:nd, d=ed(i,2); d2=d*d; l=dc([i,i+1])*d2;
     q=zeros(d,d,D); q(l(1)+1 : (d2+d+1) : l(2))=ed(i,1);
     E.data{i}=q;
  end

  E=QSpace(E);
  if ~isempty(tag), E=setitags(E,tag,3); end

end

