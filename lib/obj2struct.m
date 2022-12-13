function S=obj2struct(S)
% function S=obj2struct(S)
%
%    Suppress warning when converting object to plain structure
%
% Wb,Aug09,16

% Warning MATLAB:structOnObject :
% Calling STRUCT on an object prevents the object from hiding
% its implementation details and should thus be avoided. Use DISP
% or DISPLAY to see the visible public details of an object. 
% See 'help struct' for more information. 

  warning off MATLAB:structOnObject
  S=struct(S);
  warning on  MATLAB:structOnObject

end

