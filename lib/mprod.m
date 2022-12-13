function a=mprod(varargin)
% function a=mprod(varargin)
%
%    multiple (matrix) product
%    NB! mtimes is the matlab internal command for the operator '*'
%
% Wb,Aug11,11

  if ~nargin, a=[]; return; end
  a=varargin{1}; for k=2:nargin, a=a*varargin{k}; end

end

