function mvobj(hh,dxy)
% function mvobj(hh,dxy)
% moves given set of object handles axis by [dx dy] (units normalized)
% Wb,Feb09,10

  if nargin~=2
     eval(['help ' mfilename]);
     if nargin || nargout, wberr('invalid usage'), end, return
  end

  dxy=reshape(dxy,1,[]); n=length(dxy);

  for h=reshape(hh,1,[])
     set(h,'Units','Normalized');
   % try set(h,'Units','Normalized');
   % catch i % i=lasterror;
   %    if isequal(i.identifier,'MATLAB:class:InvalidProperty')
   %       fprintf(1,'  ERR no units property for %s objects\n',get(h,'type'));
   %    else rethrow(i); end
   % end

     p=get(h,'Position'); p(1:n)=p(1:n)+dxy;
     set(h,'Position',p);
  end

end

