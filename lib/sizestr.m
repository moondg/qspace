function s=sizestr(s,varargin)
% function s=sizestr(s [,len,sep])
% Options
%
%    sep   separator (default: 'x')
%    len   string length on the entire result (<0 is left justfied)
%
% See also num2str2(..,'--bytes')
% Wb,Feb21,21

% see also vec2str(s,'sep','x','-f');

  if isempty(s), s=''; return; end
  sep='x'; len=0;

  for i=1:nargin-1
     q=varargin{i}; if ischar(q), sep=q; else len=q; end
  end

  n=numel(s); s=matcell(reshape(s,1,[]));
  for i=1:n, s{i}=num2str(s{i}); end

  s(2,1:end-1)={sep}; s=[s{:}];

  n=[length(s), abs(len)];
  if diff(n)>0
     if len>0
          s=[repmat(' ',1,diff(n)), s];
     else s(end+1:-len)=' ';
     end
  end

end

