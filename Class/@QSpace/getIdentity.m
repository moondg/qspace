function A=getIdentity(varargin)
% function A=getIdentity([args])
%
%    same as getIdentityQS
%    but returns QSpace object (rather than plain structure).
%
% Wb,May28,16

% adapted from QSpace//contract.m

  A=class(getIdentityQS(varargin{:}),'QSpace');

end

