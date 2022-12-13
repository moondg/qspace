function [ek,Iout]=ERange(HAM,varargin)
% function [dE,Iout]=ERange(HAM [,opts])
%
%    Get energy range of terms in Hamiltonian
%    by analyzing terms in HAM.info.HH vs. distance and
%    returned as [min,max] records in ek.
%
% Options
%
%    '-s' summarize (by taking max. difference)
%
% Wb,Apr25,20

  getopt('init',varargin);
     sflag=getopt('-s');
  getopt('check_error');

  if isempty(HAM) || ~isfield(HAM.info,'HH')
     wberr('invalid usage (missing field HAM.info.HH)'); end

  HH=HAM.info.HH;
  [dk,Ik,Dk]=uniquerows([diff(HH(:,[1 3]),[],2), HH(:,[2 4])]);

  i0=find(dk(:,1)==0); n0=numel(i0);
  i1=find(dk(:,1)==1); n1=numel(i1);
  i2=find(dk(:,1)==2); n2=numel(i2);

  if n2>n1/2, i1=i2; end

  E=HAM.oez(1).op;
  A2=QSpace(getIdentityQS(E,E)); setitags(A2,{'a','b','ab'});

  ii=sort([i0;i1]); n=numel(ii); ek=zeros(n,2);
  for l=1:n, i=ii(l);
     a=HAM.ops(dk(i,2)).op;
     b=HAM.ops(dk(i,3)).op;

     H=contract(A2,'!3*',{a,'-op:a','*',{A2,b,'-op:b'}});
     ee=eigQS(H);

     xx=unique(HH(Ik{i},end)); ax=abs(xx); x=ax(find(ax==max(ax),1));

     ek(l,:)=x*[min(ee(:,1)),max(ee(:,1))];
     Iout(l)=add2struct('-',H,xx);
  end

  if sflag, ek=max(diff(ek,[],2)); end

end

