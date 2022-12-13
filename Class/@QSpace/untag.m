function A=untag(varargin)
% Usage #1: function A=untag(A)
% Usage #2: function untag(A1,A2,...)
%
%    remove itags from specified A's;
%    if no return argument is requested, the modified A's
%    are assigned back to caller by their inputname.
%
% Wb,May11,17

  for k=1:nargin, Ak=varargin{k};
     if isempty(Ak) || ~isfield(struct(Ak(1)),'info'), continue; end
     for j=1:numel(Ak)
        if isfield(Ak(j).info,'itags')
           Ak(j).info.itags=regexprep(Ak(j).info.itags,'[^\*]*','');
        end
     end
     varargin{k}=Ak;
  end

  if nargout
     if nargin~=nargout, error('Wb:ERR',...
        '\n   ERR invalid usage #1 (nargout~=nargint=%g)',nargin); end
     if nargout==1
          A=varargin{1};
     else A=varargin; end
  else
     for k=1:nargin, v=inputname(k);
        if isempty(v), error('Wb:ERR',...
           '\n   ERR invalid usage #2 (inputname not available)'); end
        assignin('caller',v,varargin{k});
     end
  end

end

