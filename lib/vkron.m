function dd=vkron(varargin)
% dd=function vkron(v1,v2,..)
%
%    kron() for vectors (similar to ikron, yet coloumns of dd
%    contain data of v1,v2, ...
%
% Wb,May11,11

% adapted from ikron.m

  if length(varargin)<2
     eval(['help ' mfilename]);
     if nargin || nargout, wbdie('invalid usage'), end, return
  end

  for i=1:nargin
     if ~isnumeric(varargin{i}) || ~isvector(varargin{i})
         wbdie('invalid usage (vector spaces required)'); end
     varargin{i}=varargin{i}(:);
  end

  dd=varargin{1};
  for i=2:nargin, d2=varargin{i}; m=size(dd,1); n=size(d2,1);
      dd=[ repmat(dd,n,1), reshape(repmat(d2,1,m)',[],1) ];
  end

end

