function [se,rr,I]=SEntropy(R,varargin)
% function [se,rr,I]=SEntropy(Rho [,eps,opts])
%
%   Calculate von Neumann entropy of the mixed (reduced)
%   density matrix R.
%
%   R may also contain non-abelian symmetries in which case
%   the output weights rr are already weighted by their 
%   corresponding multiplet degeneracies such that sum(rr)=1.
%
% Options
%
%    eps      check trace(Rho)=1 within numerical accuracy of eps (1E-12)
%   '-xm'     expand multiplet degeneracy in returned rr
%   'beta',.. assume R is a Hamiltonian => apply beta=1/T and exponentiate
%
% Wb,Jan24,11

  if nargin<1 || numel(R)~=1
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

  getopt('init',varargin);
     xmflag=getopt('-xm');
     beta=getopt('beta',0);
  eps=getopt('get_last',1E-10);

  if numel(R.Q)~=2 || ~isequal(R.Q{1},R.Q{2})
  error('Wb:ERR','invalid density matrix'); end

  if beta>0, nd=numel(R.data);
     if isdiag(R)>1 && nargout<2
        drs=get_qdim(R,2);
        for i=1:nd
           R.data{i}=exp(-beta*R.data{i});
           drs(i,2)=sum(R.data{i});
        end
        R=R/sum(prod(drs,2));
        for i=1:nd, r=R.data{i}; j=find(r);
           if ~isempty(j), rj=r(j);
              drs(i,3)=-rj*log(rj)';
           end
        end
        se=drs(:,1)'*drs(:,3);
        return
     end
     for i=1:nd, r=R.data{i}; s=size(r);
        if any(s==1),      r=exp (-beta*r);
        elseif diff(s)==0, r=expm(-beta*r);
        else error('Wb:ERR','\n   ERR invalid operator H'); 
        end
        R.data{i}=r;
     end
  end

  if isdiag(R)>1, R=diag(R); end
  if beta>0, R=R/trace(R); end

  check_trace(trace(R),eps,nargin>1);

  [rr,I]=eigQS(R);

     if size(rr,2)==2
          r1=rr(:,1).*rr(:,2);
     else r1=rr; end

     t=sum(r1); e=abs(t-1); if e>1E-8
        error('Wb:ERR','\n   ERR got tr(Rho)=%.4g !?',t); end
     i=find(r1<=0); e=norm(r1(i)); if e>1E-12
        error('Wb:ERR','invalid density matrix (%g)',e); end

  if size(rr,2)==1
     rr=sort(rr,'descend');

     i=find(rr>0);
     se=-rr(i)'*log(rr(i));

     return
  end

  if isAbelian(R)
     wblog('WRN','got cgd-dim for abelian Rho !?');
  end

  if size(rr,1)~=sum(I.DB(:,1)), error('Wb:ERR','dimension mismatch'); end

  if size(I.DB,2)>1
       dz=I.DB(:,2)./I.DB(:,1);
  else dz=ones(size(I.DB)); end

  I.rr=I.EK.data;
  I.dz=dz; n=numel(dz);

  i=find(rr(:,1)>0);
  se = rr(i,2)' * ( -rr(i,1) .* log(rr(i,1)) );

  rr=I.EK.data;

  for k=1:n
     if ~xmflag
        rr{k}=dz(k)*rr{k};
     else
        rr{k}=repmat(rr{k},1,dz(k));
     end
  end

  [rr,is]=sort(cat(2,rr{:}),'descend');

  check_trace(sum(rr),eps,nargin>1);

end

% -------------------------------------------------------------------- %

function check_trace(t,eps,fflag)

  if norm(1-t)>eps
     s=sprintf('invalid density matrix [tr(rho)=1 %+.3g]',t-1);
     if fflag || norm(t-1)>0.5, error('Wb:ERR',s);
     else 
        if norm(t-1)<1E-4, wblog(1,'WRN',s); else wblog(1,'ERR',s); end
     end
  end

end

% -------------------------------------------------------------------- %

