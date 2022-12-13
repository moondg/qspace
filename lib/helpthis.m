function helpthis
% function helpthis
% Wb,Feb23,11

  I=dbstack(1,'-completenames');
  n=numel(I);
  s={};

  for i=1:n
     if i>1 && ~isequal(I(i).file,I(i-1).file), i=i-1; break; end
     s={I(i).name, s{:}};
  end

  f=I(i).file;
  i=regexp(f,'Class/@');
  if ~isempty(i)
     f=f(i+7:end); % excluding '@' // Wb,Aug23,16
  else
     i=find(f=='/');
     if ~isempty(i), f=f(i(end)+1:end); end
  end
  if isequal(f(end-1:end),'.m'), f=f(1:end-2); end

  n=numel(s); fprintf(1,'\n');
  if n==1, help(f);
  elseif n==2 && ~isequal(f,s{end})
     n=[f '>' s{end}]; help(n);
  elseif n>1, help(f);
     if n>2, fprintf(1,'  NB! got recursive call at level=%g\n',n); end
  else wberr('invalid usage (don''t know what to do with n=%g)',n); end

end

