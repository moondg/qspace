function A=blurgauss(A,sigma)
% function A=blurgauss(A,sigma)
% 
%    Broadening of data based on normalized Gaussian kernel
%    ~ exp(-0.5*(x/sigma)^2)
%    where sigma can be chosen different for reach dimension of A;
%    absent values e.g. by having length(sigma) < ndims(A) imply that
%    sigma=0 for those dimensions, i.e. no blur. E.g, if single
%    number is specified, this blurs the first dimension only.
%
% Wb,Aug29,19

  r=ndims(A); ns=min(r,numel(sigma));

  for i=1:ns
     if sigma(i)>0, p=[i,1:i-1,i+1:r]; A= permute(A,p);
        A=blurgauss_1(A,sigma(i));     A=ipermute(A,p);
     end
  end

end

% -------------------------------------------------------------------- %
% blur first dimension

function A=blurgauss_1(A,sigma)

   s=size(A); s1=s(1); s2=prod(s(2:end));

   l=ceil(3*sigma); if l+1>s1/2, l=floor(s1/2); end
   w=exp(-((-l:l)/(sqrt(2)*sigma)).^2); w=w/sum(w); nw=numel(w);

   for k=1:s2, q=zeros(s1,nw);
      for j=1:nw, i=j-l-1;
         if i<=0
              q(1:s1+i,j) = w(j)*A(1-i:s1,  k);
         else q(i:s1,  j) = w(j)*A(1:s1-i+1,k);
         end
      end
      A(:,k)=sum(q,2);
   end

end

% -------------------------------------------------------------------- %

