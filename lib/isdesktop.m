function i=isdesktop()
% function i=isdesktop()
%
%    Check whether desktop field is set in getuser(0)->ml_env.
%
% Wb,Aug18,22

  u=get(0,'UserData');

  if ~isfield(u,'ml_env')
     if isempty(u), set_global; u=get(0,'UserData');
     else wbdie('invalid getuser(0) - empty or struct expected'); 
     end
  end

  i=u.ml_env.desktop;

end

