function displ(varargin)
% Function displ('var1','var2',val3, ...)
%
%   Display in long format
%   Arguments specified via name are looked up in calling workspace.
%
% Wb,Mar24,08

  if nargin<1
     eval(['help ' mfilename]);
     if nargin || nargout, wberr('invalid usage'), end, return
  end

  format long

  for i=1:nargin
     if ischar(varargin{i})
        evalin('caller',sprintf('disp(%s)',varargin{i}))
     else
        disp(varargin{i})
     end
  end

  format

end

