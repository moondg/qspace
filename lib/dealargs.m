function varargout=dealargs(varargin)
% Function - dealargs(varargin)
%
%    Similar to deal() but convert strings that are numbers to numbers
%    (useful when called via command line arguments that are handed
%    over as strings).
%
% Wb,Aug21,07

% if length(varargin)<1
%    eval(['help ' mfilename]);
%    if nargin || nargout, wbdie('invalid usage'), end, return
% end

  varargout=cell(1,nargout);

  for i=1:min(nargin,nargout)
      if ischar(varargin{i}), d=str2num(varargin{i}); else d=[]; end
      if isempty(d)
           varargout{i}=varargin{i};
      else varargout{i}=d; end
  end

end

