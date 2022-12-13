function varargout=cellarrsize(A,dims,varargin)
% Function cellarrsize - get size data of cell array (block matrix)
% Usage: sz = cellarrsize(A, dims, [,OPTS]);
%
% Options:
%
%   'ndim'  which dimension (default: all)
%   'all'   get dimensions for every individual submatrix 
%
% with respect to varargout - see also deal()
% Wb,Nov25,05

% mflag=1 = "matrix" = 
% cell arrary forms a matrix in the sense used by cell2mat

  if nargin>2
     getopt ('init', varargin);
        mflag = ~getopt ('all');
     if (getopt ('check_error')), return, end
  else mflag=1; end

  ndim=length(dims); sA=size(A);

  for k=1:ndim
    if mflag
       d=dims(k); n=size(A,d);
       e=[ones(1:d-1), n]; if length(e)==1, e=[e 1]; end; s=zeros(e);

       ii=prod(sA(1:d-1))*((1:n)-1)+1;

       for i=1:n, s(i)=size(A{ii(i)},d); end
    else
       d=dims(k); n=prod(size(A)); s=zeros(size(A));
       for i=1:n, s(i)=size(A{i},d); end
    end

    sz{k}=s;
  end

  varargout=sz;

return

