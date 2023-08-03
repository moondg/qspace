function q=isnumber(varargin)
% function q=isnumber(var1, var2, ...)
%
%    returns 0 if any of the input arguments
%    is not a scalar number.
%
% Wb,May14,09

  q=false(1,nargin);

  for i=1:nargin
      if isscalar(varargin{i}) && isnumeric(varargin{i})
         q(i)=true;
      end
  end

end

