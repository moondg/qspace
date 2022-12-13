function A=contract(varargin)
% function A=contract([args])
%
%    same as contractQS
%    but returns QSpace object (rather than plain structure).
%
% Wb,Sep18,06 ; Wb,Mar31,16

  A=class(contractQS(varargin{:}),'QSpace');

end

