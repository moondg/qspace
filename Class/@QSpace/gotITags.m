function [r,t] = gotITags(A)
% function r = gotITags(A)
%
%    return number of itags (if any).
%    second return argument contains listing of ilables as cell array.
%
% Wb,Aug04,12

  if nargin~=1, wblog('ERR','%s() invalid usage',mfilename); end

  r=0; t={};

  if isfield(struct(A),'info') && isfield(A.info,'itags')
     t=A.info.itags; if isempty(t), t={};  % '', [], {}
     else
        if ischar(t)
           t=strread(t,'%s','delimiter',',;');
        end
        r=numel(t);
        if r~=numel(A.Q), error('Wb:ERR',...
          '\n   ERR invalid number of itags (%d/%d)',r,numel(A.Q));
        end
     end
  end

end

