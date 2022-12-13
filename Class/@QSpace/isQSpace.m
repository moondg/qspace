function varargout=isQSpace(A,varargin)
% function [i,s]=isQSpace(A [,opts])
%
%    overloaded routine (see also MEX/isQSpace.m)
%
% see mpsIsQSpace.m for more information.
% Wb,Feb15,13

  varargout=cell(1,max(1,nargout));
  [varargout{:}]=mpsIsQSpace(struct(A),varargin{:});

end

