function i=isnumber(varargin)
% function i=isnumber(var1, var2, ...)
%
%    returns 0 if any of the input arguments
%    is not a scalar number.
%
% Wb,May14,09

  i=0;

  for k=1:nargin
      if ~isscalar(varargin{k}) || ~isnumeric(varargin{k})
      return; end
  end

  i=1;

end

