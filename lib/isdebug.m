function i=isdebug()
% function i=isdebug()
%
%    Check whether debug field is set in getuser(0)->ml_env.
%
% Wb,Nov22,23

% adapted from isbatch() // Wb,Nov22,23

% mcc-cluster jobs
  i=0; if isdeployed(), return; end

  u=get(0,'UserData');
  if isempty(u), set_global; u=get(0,'UserData'); % reload
  elseif ~isfield(u,'ml_env')
     wbwrn('invalid getuser(0) - empty or struct expected'); 
  end

  if isfield(u.ml_env,'debug') && u.ml_env.debug, i=1; end

end

