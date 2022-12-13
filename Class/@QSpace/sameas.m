function i=sameas(A,B,eps)
% function i=sameas(A,B [,eps=1E-10])
%
%    check whether QSpaces A and B are the same
%    except for minor differences such as 
%    info.cgr.cid(end) and info.cgr.nnz
%
% Wb,Aug08,16

  i=1;

  if ~isequal(size(A),size(B)), i=0; return; end
  if nargin<3, eps=1E-10; end

  for k=1:numel(A)
     if ~isequal(A(k),B(k))
        if ~isequal(A(k).Q,B(k).Q), i=0; break; end
        if ~isequal(getqdir(A(k)),getqdir(B(k))), i=0; break; end

        e=[normQS(A(k)), normQS(B(k))];
        if any(e>eps)
           e=normQS(A(k)-B(k))/max(e);
           if e>eps, i=0; break; end
        end
     end
  end

end

