function [f,I]=QSet2file(C,vflag)
% function [f,I]=QSet2file(f)
% Wb,Sep20,16

  if nargin>1
     if     isequal(vflag,'-q'), vflag=0;
     elseif isequal(vflag,'-v'), vflag=2;
     else vflag, wbdie('invalid usage'); end
  else vflag=(nargout<=1); end

  q=C.qdir; r=numel(q);
  dq=diff(q); i=find(dq); n=numel(i);
  dq=dq(i);

  if ~n && q(1)=='-'
     if vflag, wblog('WRN','using conjugate QSet'); end
     C.qdir(q=='-')='+';
  end

  if ~n || n==1 && i>=r/2 && dq==2
     f=[ C.qdir '/(' QSet2str(C,'-f;') ').cgd'];
     if nargout>1, I=[]; end
     return
  end

  ip=find(q=='+'); np=numel(ip);
  im=find(q=='-'); nm=numel(im);
  if np+nm~=r, wbdie('unexpected sorted QSet ''%s''',C.qdir); end

  I.cflag=(np<nm);
  if I.cflag
     C.qdir(q=='+')='-';
     C.qdir(q=='-')='+';
  end

  [C.qdir,I.is]=sort(C.qdir);
  q=reshape(C.qset,[],r)';

  l=find(diff(C.qdir));
  if numel(l)~=1 || size(q,1)~=r, C, q
     wbdie('unexpected rank-%d QSet',r); end

  q=q(I.is,:);
  j=1:l;   [q(j,:),i]=sortrows(q(j,:)); I.is(j)=I.is(i  );
  j=l+1:r; [q(j,:),i]=sortrows(q(j,:)); I.is(j)=I.is(i+l);

  C.qset=reshape(q',1,[]);

  if vflag && ~isequal(I.is,1:r)
     wblog('WRN','using sorted QSet'); end

  f=[ C.qdir '/(' QSet2str(C,'-f;') ').cgd'];

end

