function i=isdual_(varargin)
% function i=isdual_(Q1,Q2)
%
%    Quick and dirty check for scalar operators
%    or wave functions with all-in labels (strictly would
%    have to check for dual irreps)
%
% Wb,Oct08,19

% outsourced from blkdiag.m

  i=0;

% if ~nargin, return; end
  if nargin~=2, wbdie('invalid usage (got %g input arguments)',nargin); end

  Q1=varargin{1};
  Q2=varargin{2};

  if isequal(Q1,Q2), i=1; 
  elseif isequal(size(Q1),size(Q2))
     if isequal(sort(abs(Q1),2),sort(abs(Q2),2)), i=2; end
  end

end

