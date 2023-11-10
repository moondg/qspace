function varargout=ortho(varargin)
% function [A1,A2,Io]=ortho(Psi [,opts])
%
%    Plain wrapper routine to orthoQS()
%
% Wb,Sep05,23

  n=max(1,nargout); varargout=cell(1,n);

  [varargout{:}]=orthoQS(varargin{:});

  varargout{1}=QSpace(varargout{1}); if n>1
  varargout{2}=QSpace(varargout{2}); end

end

