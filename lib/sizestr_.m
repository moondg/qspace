function s=sizestr_(A,varargin)
% function s=sizestr_(A [,opts])
%
%    same as sizestr(), but uses size(A) of data A as input.
%
% Wb,Feb21,21

  s=sizestr(size(A), varargin{:})

end

