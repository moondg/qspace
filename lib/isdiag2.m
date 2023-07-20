function [i,ee]=isdiag2(A,eps)
% function [i,ee]=isdiag2(A[,eps])
%
%    check whether input object is diagonal matrix.
%    cell array of matrices is also accepted.
%    ee is absolute value of maximal off-diagonal element.
%
% Options
%
%   'eps'...  small (numerical) off-diagonal noise
%             still to be considered zero (0)
%
% Wb,Nov14,05 ; Wb,Jul12,11

% matlab/2015a now has its own isdiag() routine
% => renamed isdiag -> isdiag2() // Wb,Aug01,16

  if nargin<1
     helpthis, if nargin || nargout
     wbdie('invalid usage'), end, return
  end

  if iscell(A), is=zeros(size(A)); ee=zeros(size(A)); n=numel(A);
     if nargin<2
        for i=1:n
           is(i)=isequal(A{i},diag(diag(A{i})));
           if nargout>1
              ee(i)=max(reshape(abs(A{i}-diag(diag(A{i}))),[],1));
           end
           if ~is(i) && nargout<2, i=0; return; end
        end
        i=min(is);
     else
        for i=1:n
           ee(i)=max(reshape(abs(A{i}-diag(diag(A{i}))),[],1));
        end
        i=all(ee(:)<eps);
     end
     ee=max(ee(:));
  else
     if nargin<2
        i=isequal(A,diag(diag(A)));
        if nargout>1, ee=max(reshape(abs(A-diag(diag(A))),[],1)); end
     else
        ee=max(reshape(abs(A-diag(diag(A))),[],1));
        i=(ee<eps);
     end
  end

end

