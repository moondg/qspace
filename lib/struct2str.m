function s=struct2str(s,skip)
% Function: s=struct2str(s [,skip])
%
%    return given s as single line string in disp() format.
%    second argument specifies GREP pattern of lines to skip.
%
% See also param2str (which also acts on an input structure).
% Wb,Jul26,07

  if nargin==0 || ~isstruct(s)
     eval(['help ' mfilename]);
     if nargin || nargout, wbdie('invalid usage'), end, return
  end

  s=evalc('disp(s)');
  s=strread( strrep(s(1:end-2),': ','='), ...
   '%s','whitespace','\n');

  if nargin>1
     for i=length(s):-1:1
        if ~isempty(regexp(s{i},skip)), s(i)=[]; end
     end
  end

  s(:,2)={'; '}; s{end}=''; s=s';
  s=cat(2,s{:});

end

