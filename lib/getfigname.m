function s=getfigname(varargin)
% function s=getfigname([fh])
% Wb,Apr05,11

  f=[];

  while numel(varargin), v=varargin{1};
     if ishandle(v) && isequal(get(v,'type'),'figure')
        f=v; varargin(1)=[];
     else
        helpthis, if nargin || nargout
        wberr('invalid usage'), end, return
     end
  end

  if isempty(f), f=gcf; end

  for u={'fname','fout','mat'}
     s=getuser(f,u{1}); if ~isempty(s), break; end
  end

  l=length(s);
  if    ~l,   s=get(f,'tag'); if isempty(s), s='figure'; end
  elseif l>8, s=regexprep(s,'\.\w{1,3}$','');
  end

  q=regexprep(s,'[\w_.]+','');
  if ~isempty(q)
     s=regexprep(regexprep(regexprep(s,'[\[\]\{\}\(\)]+',''),':+','-'),'\s+','_');
  end

end

