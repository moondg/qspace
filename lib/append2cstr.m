function a=append2cstr(a,s)
% function a=append2cstr(a,s)
% Wb,Mar09,17

  if iscell(a)
       for i=1:numel(a), a{i}=[a{i},s]; end
  else a=[a,s];
  end

end

