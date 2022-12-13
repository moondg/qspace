function i=isvar(v)
% Function: i=isvar(v)
%
%   shortcut for exit(v,'var')
%
% Wb,Nov17,06

  i=evalin('caller', ['exist(''' v ''',''var'')']);

end

