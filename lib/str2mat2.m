function nstr = str2mat2 (str)
% nstr = str2mat2 (str)
%
%    converts a linear string with carriage return
%    into array with blank padding
%
% Wb,Jan26,04

  nstr = {};

  if    nargin<1, wblog('ERR'); eval(['help ' mfilename]); end
  if ~isstr(str), wblog('ERR'); eval(['help ' mfilename]); end

  idx = [0 find(str==sprintf('\n')) (length(str)+1)];
  n   = length(idx);

  for i=1:n-1
  nstr{i} = str(idx(i)+1:idx(i+1)-1); end

  nstr = str2mat(nstr);

return

