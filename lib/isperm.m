function q=isperm(p)
% function q=isperm(p)
%
%    Check whether input is valid permutation (1-based).
%    Return length of permutation of valid.
%
% Wb,Jul25,23

  q=0;
  if isnumeric(p), i=find(size(p)>1);
     if numel(i)<=1
        p=sort(p(:));
        if p(1)==1 && all(diff(p)==1), q=numel(p); end
     end
  end

end

