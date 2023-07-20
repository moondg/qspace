function [i,s]=isQSpace(A,varargin)
% function [i,s]=isQSpace(A [,opts])
%
%    check whether input has valid QSpace format.
%
% Options (mostly with MPS in mind)
%
%    'LRs'  specific index order
%    'sLR'  specific index order
%    'd',.. local state space dimension
%
% Wb,Aug08,08 ; Wb,Feb15,13

% NB! *this was originally MEX/mpsIsQSpace.m, and was moved to the
% QSpace class here; see also MEX/isQSpace.m which now calls *this
% Wb,Jan30,23

  sflag=(nargout>1);

% NB! by being a class routine, by construction, A is of class QSpace;
% yet, within a QSpace method, isQSpace() may be called also on
% other arguments, which then may be of type struct; this *still*
% calls *this routine iirst, i.e., before MEX/isQSpace.m
  if nargin<2
     i=1; if sflag, s=''; end
     return
  end

  i=0; s=''; 

  if ~isfield(A,'Q') || ~isfield(A,'data') || ~isfield(A,'info')
     if sflag, s=sprintf('missing fields Q, data, or info'); end
     return
  end

  n=numel(A);
  for k=1:n
     if ~iscell(A(k).Q) || ~iscell(A(k).data), if sflag
        s=sprintf('invalid {Q, data, ...} structure (%d/%d)',k,n); end
        return
     end

     if isempty(A(k).Q), q1=1; elseif isempty(A(k).Q{1}), q1=2; else q1=0; end
     if isempty(A(k).data), q2=1; else q2=0; end

     if q1 && q2, continue
     elseif xor(q1,q2), if sflag
        s=sprintf('Q or data empty, but not both (%d/%d)',k,n); end
        return
     end

     if ~isnumeric(A(k).Q{1}), if sflag
        s=sprintf('Q{1} not of type numeric (%d/%d)',k,n); end
        return
     end

     if ~isnumeric(A(k).data{1}), if sflag
        s=sprintf('data{1} not of type numeric (%d/%d)',k,n);
        return
     end
  end

  if nargin<2, i=1; return; end

  getopt('init',varargin); order=''; 
     for o={'LRs','sLR'}
        if getopt(o{1}), order=o{1}; break; end
     end

     dref=getopt('d',[]);
  getopt('check_error');

  [D,DD]=mpsGetMaxDim(A);
  id=find(D==min(D)); d=D(id(1));

  if ~isempty(dref) && ~isequal(d,dref)
     s=sprintf('local dimension mismatch (%g,%g)', d,dref);
     return
  end

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
  elseif ~isempty(order)
    wblog('WRN','invalid usage (MPS index order ''%s'')',order);
    return
  end

  i=1;

end

