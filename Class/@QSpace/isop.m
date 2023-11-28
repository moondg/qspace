function is=isop(A)
% function is=isop(A)
% Wb,Nov25,11

  is=zeros(size(A)); n=numel(A);
  for i=1:n, r=numel(A(i).Q);
     if r==2, is(i)=1;
     elseif r==3 && ...
        isfield(A(i).info,'otype') && isequal(A(i).info.otype,'operator')
        is(i)=3;
     end
  end

end

