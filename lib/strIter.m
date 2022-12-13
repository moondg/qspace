function ostr = strIter(s,n)
% Function: ss = strIter(s,n)
% 
%   repeats string s n ntimes
% 
% Wb,Nov07,05

  ostr=repmat(s,n);

% depricated

return

  if nargin==0
     eval(['help ' mfilename]);
     return
  end

  ostr={s}; ostr=ostr(ones(1,n)); ostr=cat(2,ostr{:});

return

  l=length(s); ss=''; ss(n*l)='.';

  if l>1
       for i=1:n, ss((i-1)*l+1:i*l)=s; end
  else ss(:)=s;
  end

return

