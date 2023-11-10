function [A,Iout]=appendSingletons(A,rnew,qdir)
% function A=appendSingletons(A,'new_itag[*]')
% function A=appendSingletons(A,qdir)
% function A=appendSingletons(A,r_new,qdir)
%
%    Add singleton dimension(s) to reach rank r_new.
%    qdir sets direction of new singleton dimensions (numeric; e.g. +/-1)
%
%    In usage #1 assigns the specified itag to the newly added singleton
%    dimension; this also permits to simultaneously specify qdir via a
%    trailing conjugation flag '*'.
%
%    If numel(qdir)>1, it is assumed it specifies qdir for all of A
%    where the values in qdir are ignored for the ones already existing
%    in A [this is useful if qdir was obtained via getqdir(Aother)].
%
%    NB! use makeIrop() to append irop dimension instead of this,
%    since there may be non-zero symmetry labels for all-abelians ops.
%
% Examples
%
%    A4=appendSingletons(E,'  +-')
%
% See also makeIrop, fixScalarOp.
% Wb,Sep17,18

  t0='';
  if nargin==2
     if ischar(rnew)
        if regexp(rnew,'^\s*[+-]*$')
           q=rnew; i1=find(q=='-'); i2=find(q=='+');
           rnew=numel(q);
           qdir=zeros(1,rnew); qdir(i1)=-1; qdir(i2)=1;
        else
          qdir=1;
          t0=regexprep(rnew,'\*+$(?@qdir=-1;)','');
        end
     else qdir=rnew;
     end
     if numel(qdir)==1, qdir=[getqdir(A), qdir]; end
     rnew=numel(qdir);

  elseif nargin==3
     n=numel(qdir);
     if n==1, qdir=repmat(qdir,1,rnew);
     elseif numel(qdir)~=rnew
        wbdie('invalid qdir (len=%g/%g)',numel(qdir),rnew);
     end
  end

  if nargout>1
     Iout=add2struct('-',rnew,qdir);
  end

  for k=1:numel(A)
     Ak=A(k); qda=getqdir(Ak);
     r=numel(Ak.Q); if r==rnew, continue; end
     if r>rnew, wbdie('invalid rank r=%g->%g !?',r,rnew); end
     E=getvac(Ak);

     for l=r+1:rnew
        tl=regexprep(getitags(Ak,l-1),'\*+$','');
        A2=QSpace(getIdentityQS(Ak,l-1,E,[1 3 2]));
        A2=setitags(A2,{tl,tl,t0});

        if qdir(l)>0
             if qda(l-1)>0, ic=2; else ic=1; end
        else if qda(l-1)>0, ic=1; else ic=2; end; ic=[num2str(ic) '*'];
        end
        Ak=contract(Ak,l-1,A2,ic);
        qda(l)=qdir(l);
     end
     A(k)=Ak;
  end

end

