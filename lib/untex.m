function s=untex(s)
% function s=untex(s)
% Wb,Nov11,09

% s=strrep(s,'_','\_');
  s=regexprep(s,'([^\\])_','$1\\_');

end

