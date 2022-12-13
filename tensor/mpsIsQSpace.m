function [i,s]=mpsIsQSpace(A,varargin)
% function [i,s]=mpsIsQSpace(A [,opts])
%
%    check whether input A is a valid MPS in QSpace format.
%    also used by isQSpace.m (which itself is overloaded by @QSpace).
%
% Options
%
%    'LRs'  specific index order
%    'sLR'  specific index order
%    'd',.. local state space dimension
%
% Wb,Aug08,08 ; Wb,Feb15,13

  if ~nargin
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  order=[]; i=0; s='';
  if isa(A,'QSpace') && nargin<2 && nargout<2, i=1; return, end

  getopt ('init', varargin);
      o='LRs'; if getopt(o), order=o; end
      o='sLR'; if getopt(o), order=o; end
      dref=getopt('d',[]);
  getopt('check_error');

  if ~isa(A,'QSpace') && (~isfield(A,'Q') || ~isfield(A,'data'))
     s=sprintf('not of type {Q, data, ...} structure'); return
  elseif isempty(A), s=sprintf('got empty object'); return
  end

  n=numel(A);
  for k=1:n
     if (isempty(A(k).Q) || isempty(A(k).Q{1})) && isempty(A(k).data)
     i=1; return; end

     if xor(isempty(A(k).Q), isempty(A(k).data))
        s=sprintf('Q or data empty, but not both (%d/%d)', k,n);
        return
     end

     if ~iscell(A(k).Q) || ~iscell(A(k).data)
        s=sprintf('invalid {Q, data, ...} structure (%d/%d)', k,n);
        return
     end

     if ~isnumeric(A(k).Q{1})
        s=sprintf('Q{1} not of type numeric (%d/%d)', k,n);
        return
     end

     if ~isnumeric(A(k).data{1})
        s=sprintf('data{1} not of type numeric (%d/%d)', k,n);
        return
     end
  end

if nargin<2, i=1; return; end

  [D,DD]=mpsGetMaxDim(A);
  id=find(D==min(D)); d=D(id(1));

  if ~isempty(dref) && ~isequal(d,dref)
     s=sprintf('local dimension mismatch (%g,%g)', d,dref);
     return
  end

  if ~isempty(order)

     if isequal(order,'LRs')

        if all(id~=3) || any(DD(:,3)>d)
           s=sprintf('not of LRs order (id=%s, d=%d ?)', ...
           vec2str(id), max(DD(:,3))); return
        elseif DD(1,1)~=1
           s=sprintf('not of LRs order (dl=%d)', DD(1,1));
           return
        elseif DD(end,2)~=1
           s=sprintf('not of LRs order (dr=%d)', DD(end,2));
           return
        end

     elseif isequal(order,'sLR')

        if all(id~=1) || any(DD(:,1)>d)
           s=sprintf('not of LRs order (id=%s, d=%d ?)',...
           vec2str(id), max(DD(:,1))); return
        elseif DD(1,2)~=1
           s=sprintf('not of LRs order (dl=%d)', DD(1,2));
           return
        elseif DD(end,3)~=1
           s=sprintf('not of LRs order (dr=%d)', DD(end,3));
           return
        end

     else
        wblog('ERR','Invalid argument #2 (LRs order)');
        order, return
     end
  end

  i=1;

end

