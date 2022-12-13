function s=hostname
% Usage: n=hostname()

  s=evalc('!hostname');
  s=strtrim(s);

end

