function s=vsprintf(x,sep,fmt)
% function s=vsprintf(x,sep,fmt)
%
%    Print vector with specified print format and separators
%
% see also vec2str.m
% Wb,Feb07,18

  if nargin<2, sep=' ';  end
  if nargin<3, fmt='%g'; end

  n=numel(x); s=repmat({''},2,n);
  for i=1:n, s{1,i}=sprintf(fmt,x(i)); end
  s(2,1:end-1)={sep};

  s=[s{:}];

end

