function s=structdisp(varargin)
% Function: s=structdisp(s1, s2, ... [,OPTS])
%
%    Aligned display of set of structures
%    based on MatLabs native display of structures
%
% Options
%
%   '-x',..   GREP pattern of lines to skip
%   '-h'      horizontal display
%   '-dn',..  number of columns to indent whole display
%             (vertical mode only)
%
% Wb,Dec20,07

  if ~nargin || ischar(varargin{1})
     if ~helpthis(nargout,varargin{:}), wbdie('invalid usage'); end
     return
  end

  vflag=1;

  getopt('init',varargin);
     xpat =getopt('-x',[]);
     hflag=getopt('-h');
     if getopt('-q'), vflag=0; end
     dn   =getopt('-dn',[]);
     istr =getopt('istr',[]);
  varargin=getopt('get_remaining');

  nargs=length(varargin); s={}; mark=zeros(1,nargs);
  for i=1:nargs
     if     isstruct(varargin{i}), s{end+1}=disp_1(varargin{i});
     elseif isobject(varargin{i}), s{end+1}=disp_1(struct(varargin{i}));
     else
        wbdie('invalid usage (not of type struct)');
     end
  end

  s(2,1:end-1)={{'--','--'}};
  s=cat(1,s{:}); n=size(s,1);

  if ~isempty(xpat), mark=zeros(1,n);
     for i=1:n
        if ~isempty(regexp(s{i,1},xpat)), mark(i)=1; end
     end
     s(find(mark),:)=[]; n=size(s,1);
  end

  if hflag
     s(:,3)={'='}; s=s(:,[1 3 2]);
     s(:,4)={', '}; s{end,4}=''; s=s';

     s=regexprep([s{:}],'[,;]*[ ]*--[ :=]*--[,;]*[ ]*',';  ');
     if ~isempty(istr), s=[istr ': ' s]; end
  else
     s(:,3)={': '}; s=s(:,[1 3 2]);

     nn=zeros(n,1);
     for i=1:n, nn(i)=length(s{i}); end
     if isempty(dn), dn=max(6,max(nn)+3); end

     for i=1:n
        s{i}=[blanks(dn-nn(i)), s{i,:}];
     end

     s=[' '; s(:,1); ' ']';
     if ~nargout && ~isempty(istr), s={' ', [ '   ' istr ' ='], s{:}}; end

     s(2,:)={ char(10) }; s{2,end}='';

     s=regexprep([s{:}],'--[ :=]*--','');
  end

  if ~nargout
     if nargin && vflag, v=inputname(1);
        if ~isempty(v)
           v={v,' ='}; if nargin>1, v{2}=':'; end
           fprintf(1,'\n   %s%s',v{:});
        end
     end
     disp(s); clear s
  end

end

% -------------------------------------------------------------------- %
% structdisp(S) adapts the output of matlab's native disp(S)
% Wb,Nov28,23

function s=disp_1(S)
  s=evalc('disp(S)');

  s=regexprep(s(1:end-1),': ','\n');
  s=reshape(strread(s,'%s','whitespace','\n'),2,[])';

  for i=1:size(s,1)
     if ~isempty(regexp(s{i,2},'^.\d+-D '))
        sz=sprintf('x%d',size(getfield(S,s{i,1})));
        s{i,2}=sprintf('%-20s%s',s{i,2},sz(2:end));
     end
  end

end

% -------------------------------------------------------------------- %
