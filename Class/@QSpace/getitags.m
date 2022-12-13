function t=getitags(A,it)
% function t=getitags(A)
%
%    return itags as cell array.
%
% adapted from QSpace/gotITags.m
% Wb,Apr10,14

  if nargin<1 || nargin>2 || numel(A)~=1
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  t={};

  if isfield(struct(A),'info') && isfield(A.info,'itags')
     t=A.info.itags; if isempty(t), t={};  % '', [], {}
     else
        if ischar(t)
           t=strread(t,'%s','delimiter',',;')'; % ';*;*'
        end
        r=numel(t);
        if r~=numel(A.Q), error('Wb:ERR',...
          '\n   ERR invalid number of itags (%d/%d)',r,numel(A.Q));
        end
     end
  end

  if nargin>1
     if numel(it)==1, t=t{it}; else t=t(it); end
  end

end

