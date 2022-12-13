function q=setbatch(bval)
% function q=setbatch([val])
%
%    Set batchmode to set(0,'UserData').
%    To unset batch mode, use setbatch(0).
%    Return value if requested contains old value.
%
% Wb,Sep06,07

  if ~nargin, bval=1;
  elseif ischar(bval)
     q=str2num(bval); if ~isempty(q), bval=q; end
  end

  u=get(0,'UserData');
  if isempty(u), set_global; u=get(0,'UserData');
  elseif ~isfield(u,'ml_env')
     wbdie('invalid getuser(0).ml_env'); disp(u);
  end

  if nargout, q=u.ml_env.batch; end

  u.ml_env.batch=bval;
  if bval>0,    u.ml_env.desktop=0;
  elseif ~bval, u.ml_env.desktop=set_global('--test-dkt');;
  end

  set(0,'UserData',u);

  if ~nargout, clear q; end

end

