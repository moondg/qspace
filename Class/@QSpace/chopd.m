function A=chopd(A,aref)
% function A=chopd(A [,aref])
%
%    apply chopd() to data
%    => sets small matrix elements to zero
%    => removes imaginary noise if any
%
% Wb,Oct17,11

  a2=zeros(size(A));
  for k=1:numel(A), a2(k)=normQS(A(k)); end

  if nargin<2, aref=max(a2(:)); end
  eps=1E-12*aref;
  aps=1E+02*aref;

  for k=1:numel(A)
     Ak=A(k); m=length(Ak.data); if m==0, continue; end
     mark=zeros(m,1);

     for i=1:m
        if norm(Ak.data{i}(:))<eps, mark(i)=1;
        else
           Ak.data{i}=((Ak.data{i}+aps)-aps);
           if ~isreal(Ak.data{i})
              if norm(imag(Ak.data{i}(:)))<eps
              Ak.data{i}=real(Ak.data{i}); end
           end
        end
     end

     I=find(mark); if ~isempty(I)
        for i=1:length(Ak.Q), Ak.Q{i}(I,:)=[]; end
        Ak.data(I)=[];
     end

     A(k)=Ak;
  end

end

