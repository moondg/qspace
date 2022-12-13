function s=hostid()
% function s=hostid()
% Wb,Apr19,14

% s=system_getline('~/bin/hostid.pl');
  s=evalc('!hostname'); s=strtrim(s);

end

