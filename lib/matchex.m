function i=matchex(varargin)
% function i=matchex(['-i',] s1,regexp [,opts])
%
%    check match of string using regexp(i)
%
% Options
%
%   '-i'  ignore case
%
% Wb,Jul25,18

% renamed from matches -> matchex to avoid startup warning
% with matlab/2020b which also introduced simular function
% Wb,Oct12,20

  if nargin>2
     if ~isequal(varargin(1:end-2),{'-i'})
        wberr('invalid usage'); end
     i=~isempty(regexpi(varargin{2:end}));
  else
     i=~isempty(regexp (varargin{:}));
  end

end

