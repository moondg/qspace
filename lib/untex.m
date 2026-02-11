function s=untex(s,varargin)
% function s=untex(s [,regexprep1, regexrep2, ...])
% Wb,Nov11,09 ; Wb,2025

  s=regexprep(s,'([^\\])_','$1\\_');

  for i=2:nargin, q=varargin{i-1}; i=[];
     if ischar(q), q=textscan(q,'%s','whitespace','/'); q=q{1}; end
     if iscell(q) && ~isempty(q) && isequal(q{1},'s'), q(1)=[]; end
     if ~iscell(q) || numel(q)<2, wbdie('invalid usage'); end
     s=regexprep(s,q{1:2});
  end

end

