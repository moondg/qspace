function s=param2str(param,varargin)
% Function: param2str - convert parameter structure to string
% Usage: s=param2str(param [, OPTS, replace strings])
%
%    param is the parameter structure
%    replace strings are given by ..., 'str1', 'str1r', 'str2', 'str2r', ...
%    where str1 would be replace by str1r and so on.
%
% Options
%
%   '--tex'   escale some Greek letters and underscore
%   '-x',..   exclude fields matching grep pattern
%   'sep',..  field separator (', ')
%   'fmt',..  format for numbers and vectors
%
% See also struct2str.m
% Wb Jan 2007

  getopt ('init', varargin);
     sep = getopt('sep',', ');
     fmt = getopt('fmt',  {});

     tex = getopt('--tex'   );
     xpat= getopt('-x',   {}); if ~iscell(xpat), xpat={xpat}; end

     nrgflag=getopt('--nrg');
     if nrgflag, xpat{end+1}='nrgIO|ALambda|wsys'; end

  varargin=getopt('get_remaining');

  if isempty(param), s=''; return; end

  fn=fieldnames(param)';
  istr='';

  i=find(cellfun(@(x) ~isempty(x), regexpi(fn,'istr|info')));
  if numel(i)==1
     if ~isempty(xpat) && isempty(regexp(fn{i},xpat))
        istr=getfield(param,fn{i});
     end
     fn(i)=[];
  end
  if ~isempty(xpat)
     i=find(cellfun(@(x) ~isempty(x), regexp(fn,strjoin(xpat,'|'))));
     fn(i)=[];
  end

  m=length(fn);
  fn(4,:)={sep}; fn(2,:)={'='};

  if ~isempty(fmt), fmt={'fmt',fmt}; end

  for i=1:m
     d=getfield(param,fn{1,i});
   % if iscell(d) && numel(d)==2 && isequal(d{1},'--const')
   %  % e.g., oJ2: {'--const', 0.5000} // Wb,Dec22,21
   %    fn{1,i}=regexprep(fn{1,i},'^o([A-Z])','$1');
   %    d=d{2};
   % end
     if ischar(d)
          fn{3,i}=[ '''' d '''' ];
     elseif ~isnumeric(d), fn(:,i)={''}; continue
     elseif numel(d)>4
          fn{3,i}=sprintf('[%dx%d array]',size(d));
     elseif numel(d)>1
          fn{3,i}=[ '[' vec2str(d,fmt{:}), ']'];
     elseif ~isempty(fmt)
          fn{3,i}=vec2str(d,fmt{:});
     else fn{3,i}=sprintf('%.4g',d);
     end
  end

  for i=m:-1:1
     fn{4,i}='';
     if ~isempty(fn{1,i}), break; end
  end

  s=[fn{:}];
  if istr, s=[istr ': ' s]; end

  if tex, s=str2tex(s);
     if nrgflag
        symstr = @(x) format_sym_str(x);

        s=regexprep(s,'\<Etrunc\>','E_{tr}');
        s=regexprep(s,'\<(1|one)[\s-]*channel[\s-][^,]*\((SIAM)\)','$2','ignorecase');
        s=regexprep(s,'sym=''([^'']+)''','${symstr($1)}');
     end
  end

  for i=2:2:length(varargin)
     s=strrep(s,varargin{i-1},varargin{i});
  end

end

function s=format_sym_str(s)
   s=regexprep(s,'\<SU2spin\>','SU2_{sp}');
   s=regexprep(s,'\<SU2charge\>','SU2_{ch}');
   s=regexprep(s,',','{\\otimes}');
end

