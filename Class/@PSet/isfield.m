function i=isfield(P,f)
% Function i=isfield(P,f)
%
%    check whether field f is part of *this PSet.
%    returns field-id if field exists, 0 otherwise.
%
% Wb,Aug26,08 ; Wb,Nov22,11

% tags: idfield, fieldid

% NB! use builtin('isfield',...) for this purpose within class member
% routines (as referring to specific fields is not of interest outside
% the class object) -- Wb,Nov22,11
% MatLab 2008a: isfield(P,'field') always returns 0
% MatLab 2008a no longer allows to access object fields this way !??
% Wb,Aug26,08
% i=isfield(struct(P),f); % isfield(P,f);

  n=numel(P.vars); i=0;
  for j=1:n, if isequal(P.vars{j},f), i=j; break; end, end

end

