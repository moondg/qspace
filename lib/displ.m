function displ(varargin)
% Function displ('var1','var2',val3, ...)
%
%   Display in long format
%   Arguments specified via name are looked up in calling workspace.
%
% Wb,Mar24,08

  if nargin<1
     eval(['help ' mfilename]);
     if nargin || nargout, wbdie('invalid usage'), end, return
  end

  format longg

  for i=1:nargin, x=varargin{i};
     if ischar(x)
        evalin('caller',sprintf('disp(%s)',x))
     elseif numel(x)==1 && isfloat(x)
        disp(sprintf('\n   %.16g\n',x));
     else
        disp(x)
     end
  end

  format

end

