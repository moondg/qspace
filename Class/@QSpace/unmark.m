function varargout=unmark(varargin)
% Usage #1: function A=unmark(A)
%
%    Remove all marks from itags, where marks are considered 
%    trailing special characters other than conjugate flag `*'
%    such as [~^'+-].
%
% Usage #2: function unmark(A1,A2,...)
%
%    With no return argument requested, all modified A's
%    are assigned back to caller by their inputname
%    (hence, in this case input arguments cannot be r-values
%    out of exressions, but must be actual variables).
%
% See also QSpace/untag.m
% Wb,Jun29,23

% adpated from QSpace/untag.m // Wb,Jun29,23
  tpat='[~^''?+-]*(\*?)$';

  if nargout
     if nargin~=nargout, wbdie(...
        'invalid usage (nargout / nargin = %d / %d)',nargout,nargin); end
     varargout=varargin;
  end

  for k=1:nargin, Ak=varargin{k}; m=0;
     if isempty(Ak) || ~isfield(struct(Ak(1)),'info'), continue; end
     for j=1:numel(Ak)
        if isfield(Ak(j).info,'itags')
           t0=Ak(j).info.itags; t2=regexprep(t0,tpat,'$1');
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

