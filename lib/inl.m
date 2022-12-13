function inl(n)
% Function : inl - inserts emty lines on command prompt
% Usage    : inl nr or inl(nr)
%
% Wb,Mar25,05

  if nargin~=1, eval(['help ' mfilename]); return; end
  if isstr(n), n=str2num(n); end

  fprintf(1,repmat(char(10),1,n));

end

