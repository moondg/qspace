function setfig(fs)
% function setfig(fs)
% Wb,Apr04,13

% adapted from setax()

  if nargin~=1 || nargout
     helpthis, if nargin || nargout, wberr('invalid usage'), end
     return
  end

  if ~ischar(fs)
     if numel(fs)~=1 || (~isnumeric(fs) && ~isfig(fs))
        wberr('invalid figure handle');
     end
  else
     h=findall(groot,'type','figure','tag',fs);
     if numel(h)~=1
        if numel(h)>1
             wberr('%g figures found with tag ''%s''',fs);
        else wberr('no figure with tag ''%s'' yet',fs); 
        end
     end
     fs=h;
  end

  set(0,'CurrentFigure',fs);

end

