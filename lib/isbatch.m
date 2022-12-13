function i=isbatch()
% function i=isbatch()
%
%    Check whether batch field is set in getuser(0)->ml_env.
%
% Wb,Sep06,07 ; Wb,Aug18,22

% mcc-cluster jobs // Wb,Feb15,17
  i=isdeployed(); if i, return; end

  u=get(0,'UserData');
  if isempty(u), set_global; u=get(0,'UserData');
  elseif ~isfield(u,'ml_env')
     wbdie('invalid getuser(0) - empty or struct expected'); 
  end

  i=u.ml_env.batch;

end

