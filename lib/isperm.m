function q=isperm(p)
% function q=isperm(p)
%
%    Check whether input is valid permutation (1-based).
%    Return length of permutation of valid.
%
% Wb,Jul25,23

  q=0;
  if isnumeric(p), i=find(size(p)>1);
     if numel(i)<=1 % enforce vector
        p=reshape(p,1,[]); r=numel(p);
        if isequal(p,round(p)) && all(p>=1 & p<=r) && all(diff(sort(p)))>=1
           q=r;
        end
     end
  end

end

