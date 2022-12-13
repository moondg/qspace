function label3(xstr,ystr,zstr,tstr,zorient)
% Function label3(xlabel,ylabel,zlabel,tstr [,zorientation])
%
%    adds x/y/z label and title to current figure
%
% Wb,Feb28,01

  if nargin>0, xlabel(xstr); end
  if nargin>1, ylabel(ystr); end
  if nargin>2, zlabel(zstr); end
  if nargin>3, title (tstr); end

  if nargin>4
     zorient = num2str(zorient);
     if strcmp(zorient,'y') || strcmp(zorient,'Y') || strcmp(zorient,'1')
          set(get(gca, 'zlabel'), 'Rotation', 0)
     else set(get(gca, 'zlabel'), 'Rotation', 90)
     end
  end

end

