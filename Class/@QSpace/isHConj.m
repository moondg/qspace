function q=isHConj(A,deps)
% function q=isHConj(A [,deps])
%
%    Check whether input A represents Hermitian conjugate operator
%    within treshold deps (1E-12).
% 
% Wb,Apr04,23

  if nargin<2, deps=1E-12; end
  q=false(size(A));
  for k=1:numel(A), Q=A(k).Q; rk=numel(Q);
     if ~rk, nd=numel(A(k).data);
        if ~nd, q(k)=1;
        elseif nd==1
           x=A(k).data{1}; s=size(x);
           if numel(s)>2 || diff(s), continue; end
           e=norm(x-x'); if e<=deps, q(k)=1; end
        end
     elseif rk==2 || rk==3 && ~norm(Q{3})
        depsk=max(1,normQS(A(k)))*eps;
        nk=numel(A(k).data); qk=zeros(1,nk);
        for i=1:nk, x=A(k).data{i}; s=size(x);
           if numel(s)>2 || diff(s), continue; end
           e=norm(x-x'); if e<=depsk, qk(i)=1; end
        end
        if all(qk), q(k)=1; end
     end
  end

end

