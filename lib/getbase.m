function varargout=getbase(varargin)
% Function [a1,a2,...]=getbase('var1','var2',...)
% Wb,Sep27,07

% WARNING, Wb,Feb14,11 ----------------------------------------------- %
% when called in batch-mode which requires function
% base may return empty for most parts! // tags: ISBATCH
% wblog('WRN','calling getbase(%s) in batch-mode', strhcat(varargin));
% -------------------------------------------------------------------- %

% NB! calling caller from caller acts cyclic between this and caller
% evalin('caller','evalin(''caller'',''whos, a=3;'')');
% would set a=3 in *THIS* workspace (caller->caller->caller (see i ml)

  global varx__
  evalin('base','global varx__');

  if nargout
     varargout=cell(1,max(nargout,nargin));
     for i=1:nargin, v=varargin{i};
        evalin('base',sprintf(...
        'if exist(''%s'',''var''), varx__=%s; else varx__=[]; end',v,v))
        varargout{i}=varx__;
     end
  else
     for i=1:nargin, v=varargin{i};
        evalin('base',sprintf(...
           'if exist(''%s'',''var''), varx__=%s; else varx__=[]; end',v,v))
        if ~isbatch
           assignin('caller',v,varx__);
        else
           x=varx__; varx__=0; if isempty(x)
              evalin('caller',sprintf(...
              'if exist(''%s'',''var''), global varx__; varx__=1; end',v))
           end
           if ~isempty(x) || ~varx__, assignin('caller',v,x); end
        end
     end
  end

  clear global varx__

end

