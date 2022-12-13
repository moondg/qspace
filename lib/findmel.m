function [i, j, dd] = findmel (M, mij, delta)
% Function : findmel - find matrix elements
% Usage    :
%
%     [i, j, dd] = findmel (M, mij [, delta])
%     (delta<0 finds closest value(s)!)
%
% Wb,Oct03,02, May10,03 - using [i j dd] = find (<condition>)

  if nargin==0 | nargin>4
     help findmel
     return
  end

  if ~exist('mij'); mij = 0; end
  if isempty(M), i=[]; j=[]; dd=[]; return; end;

  if ~exist('delta')
      [i j dd] = find(M==mij);
  else
      if delta>0
           [i j dd] = find( abs(M-mij) < delta );
      else
        if delta==-1
           [i j dd] = find( abs(M-mij) == min(min(abs(M-mij))) );
        else
           [xx,ii]=sort(abs(M(:)-mij)); nel=min([-delta,length(xx)]);
           [i j dd] = find( abs(M-mij) <= xx(nel) );
        end
      end
  end

  [n1 n2] = size(M);
  dd=M(i+(j-1)*n1);

  if isvector(M)
     if size(M,2)>1, i=j; end
     j=dd; clear dd;
  end

return

