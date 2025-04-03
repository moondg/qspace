function s=param2str(param,varargin)
% function s=param2str(param [,opts, replace strings])
%
%    Convert parameter structure <param> to string.
%    Replace strings are given by ..., 'str1', 'str1r', 'str2', 'str2r', ...
%    where str1 would be replace by str1r and so on.
%
% Options
%
%   '--tex'   make use of Greek letters, underscores, and exponentials
%   '-x',..   exclude fields matching grep pattern or values like [], 0, etc.
%             this also accepts cell input, like '-x',{some_regexp, [], 0}.
%   '-m',..   only include fields that match given grep pattern
%   'sep',..  field separator (', ')
%   'fmt',..  format for numbers and vectors
%
% See also struct2str.m
% Wb Jan 2007

% Examples // see also mat_to_headerSE.m
% header('SW',param2str(param,'--nrg','-x','H0str','--tex'),{'FontSize',10});

  getopt ('init', varargin);
     sep = getopt('sep',', ');
     fmt = getopt('fmt',  {});

     tex = getopt('--tex');
     fpat= getopt('-m',{}); if ~iscell(fpat), fpat={fpat}; end
     xpat= getopt('-x',{}); if ~iscell(xpat), xpat={xpat}; end

     nrgflag=getopt('--nrg');
     if nrgflag, xpat{end+1}='nrgIO|ALambda|wsys'; end

  varargin=getopt('get_remaining');

  if isempty(param), s=''; return; end

  fn=fieldnames(param)';
  istr='';

  if ~isempty(xpat), n=numel(xpat); ix=[];
     for l=1:n
        if isempty(xpat{l}) || ~ischar(xpat{l}), ix(end+1)=l;
        end
     end
     if ~isempty(ix), nf=numel(fn); q=zeros(1,nf);
        ff=cell(1,nf); for i=1:nf, ff{i}=getfield(param,fn{i}); end
        for l=ix
           for i=1:nf
              if ~q(i) && isequal(ff{i},xpat{l}), q(i)=l; end
           end
        end
        xpat(ix)=[];
        i=find(q); if ~isempty(i), fn(i)=[]; clear ff; end
     end
  end

  i=find(~cellfun(@isempty, regexpi(fn,'istr|info')));
  if numel(i)==1
     if ~isempty(xpat) && isempty(regexp(fn{i},xpat))
        istr=getfield(param,fn{i});
     end
     fn(i)=[];
  end
  if ~isempty(xpat)
     i=find(~cellfun(@isempty, regexp(fn,strjoin(xpat,'|'))));
     fn(i)=[];
  end
  if ~isempty(fpat)
     i=find( cellfun(@isempty, regexp(fn,strjoin(fpat,'|'))));
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
     s=regexprep(s,...
       '(\d[\.\d])*[eE]([+-]\d+)',['$1' char(1) '10^{' char(2) '$2}']);
     if ~isempty(s==1)
        s=regexprep(s,'([^0-9\.])1\x01','$1');
        s=regexprep(s,'\x01','{\\cdot}');
        s=regexprep(s,'\x02\+0*','' );
        s=regexprep(s,'\x02\-0*','-');
     end

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

