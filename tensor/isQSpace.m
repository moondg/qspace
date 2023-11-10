function [i,s]=isQSpace(A,varargin)
% function [i,s]=isQSpace(A [,opts])
%
%   wrapper routine that calls QSpace/isQSpace.m
%
% formerly called mpsIsQSpace.m
% Wb,Jan30,23

  if ~nargin
     helpthis, if nargin || nargout
     wbdie('invalid usage'), end, return
  end

  order=''; i=0; s='';

  if isa(A,'QSpace')
     if nargin<2 && nargout<2, i=1; return, end
  elseif ~isfield(A,'Q') || ~isfield(A,'data') || ~isfield(A,'info')
     if nargout>1, s=sprintf('missing fields Q, data, or info'); end
     return
  elseif isempty(A)
     if nargout>1, s=sprintf('empty object'); end
     return
  end

  varargout=cell(1,max(1,nargout));
  [varargout{:}]=isQSpace(QSpace(A),varargin{:});

end

