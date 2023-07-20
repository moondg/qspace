function varargout=untag(varargin)
% Usage #1: function A=untag(A)
%
%    Remove itags from specified QSpace A
%    i.e., all characters other than trailing conjugate flag `*'.
%
% Usage #2: function untag(A1,A2,...)
%
%    With no return argument requested, all modified A's
%    are assigned back to caller by their inputname
%    (hence, in this case input arguments cannot be r-values
%    out of exressions, but must be actual variables).
%
% See also QSpace/unmark.m
% Wb,May11,17

  if nargout
     if nargin~=nargout, wbdie(...
       'invalid usage (nargout / nargin = %d/%d)',nargout,nargin); end
     varargout=varargin;
  end

  tpat='[^\*]*';

  for k=1:nargin, Ak=varargin{k}; m=0;
     if isempty(Ak) || ~isfield(struct(Ak(1)),'info'), continue; end
     for j=1:numel(Ak)
        if isfield(Ak(j).info,'itags')
           t0=Ak(j).info.itags; t2=regexprep(t0,tpat,'');
           if ~isequal(t0,t2)
              Ak(j).info.itags=t2; m=m+1;
           end
        end
     end
     if ~m, continue
     elseif nargout, varargout{k}=Ak;
     else
        v=inputname(k); if isempty(v), wbdie(['invalid usage #2 ' ... 
         '(inputname for argument %d not available)'],k); end
        assignin('caller',v,Ak);
     end
  end

end

