function q=isint(x,deps)
% function q=isint(x [,deps])
% Wb,Oct11,21

  q=false;

  if ~isnumeric(x) || isempty(q)
  elseif isinteger(x), q=true;
  elseif numel(x)==1 && nargin<2, q=~(x-round(x));
  else
     x=x(:); e=norm(x-round(x));
     if ~e, q=true;
     elseif max(abs(x))>=0.1
        if nargin<2, deps=1E-14; end
        if e<deps, q=true; end
     end
  end

end

