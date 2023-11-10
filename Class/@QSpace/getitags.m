function [t,cflag,mark]=getitags(A,it)
% function [t,conj,mark]=getitags(A [,it])
%
%    return itags as cell array, except if single itag is requested
%    via the index set it = index of tags. In the later case,
%    if two additional arguments are requested, trailing conj-flag (*)
%    and marker (') are returned as counts in cflag and mark,
%    respectively.
%
% adapted from QSpace/gotITags.m
% Wb,Apr10,14

  if nargin<1 || nargin>2 || numel(A)~=1
     helpthis, if nargin || nargout
     wbdie('invalid usage'), end, return
  end

  t={};

  if isfield(struct(A),'info') && isfield(A.info,'itags')
     t=A.info.itags; if isempty(t), t={};  % '', [], {}
     else
        if ischar(t)
           t=strread(t,'%s','delimiter',',;')'; % ';*;*'
        end
        r=numel(t);
        if r~=numel(A.Q)
           wbdie('invalid number of itags (%d/%d)',r,numel(A.Q));
        end
     end
  end

  if nargin>1
     if numel(it)~=1, t=t(it);
     else t=t{it};
        if nargout>1
           n=0; t=regexprep(t,'(\*+)$(?@n=numel($1);)','');
           cflag=mod(n,2);
		if nargout>2
		   n=0; t=regexprep(t,'(''+)$(?@n=numel($1);)','');
		   mark=mod(n,2);
		end; end
     end
  end

end

