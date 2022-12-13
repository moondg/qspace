function s=structdisp(varargin)
% Function: s=structdisp(s1, s2, ... [,OPTS])
%
%    aligned display of set of structures
%    based on MatLabs native display of structures
%
% Options
%
%   '-x',..   GREP pattern of lines to skip
%   '-h'      horizontal display
%   '-dn',..   number of columns to indent whole display
%             (vertical mode only)
%
% Wb,Dec20,07

  if ~nargin
     eval(['help ' mfilename]);
     if nargin || nargout, wberr('invalid usage'), end, return
  end

  getopt('init',varargin);
     xpat =getopt('-x',[]);
     hflag=getopt('-h');
     dn   =getopt('-dn',[]);
     istr =getopt('istr',[]);
  varargin=getopt('get_remaining');

  narg=length(varargin); s=cell(1,narg); mark=zeros(1,narg);
  for i=1:narg
     if ~isstruct(varargin{i}) && ~isobject(varargin{i})
        wberr('invalid usage (not of type struct)'); end

     if ~isempty(varargin{i}) && ~isempty(fieldnames(varargin{i}))
        s{i}=sprintf('disp(varargin{%g}); ',i);
     else mark(i)=1; end
  end

  s(find(mark))=[];
  s(2,:)={'disp(''-- : --''); '}; s{2,end}='';

  s=evalc(cat(2,s{:})); s=[char(10) s(1:end-1)];
  s=reshape(strread(regexprep(s,': ','\n'),'%s','whitespace','\n'),2,[])';

  if ~isempty(xpat), n=size(s,1); mark=zeros(1,n);
    for i=1:n
      if ~isempty(regexp(s{i,1},xpat)), mark(i)=1; end
    end
    s(find(mark),:)=[];
  end

  if hflag
     s(:,3)={'='}; s=s(:,[1 3 2]);
     s(:,4)={', '}; s{end,4}=''; s=s';

     s=regexprep(cat(2,s{:}),'[,;]*[ ]*--[ :=]*--[,;]*[ ]*',';  ');

     if ~isempty(istr), s=[istr ': ' s]; end
  else
     s(:,3)={': '}; s=s(:,[1 3 2]);

     n=size(s,1); nn=zeros(n,1);
     for i=1:n, nn(i)=length(s{i}); end
     if isempty(dn), dn=max(6,max(nn)+3); end

     for i=1:n
        s{i}=[blanks(dn-nn(i)), s{i,:}];
     end

     s=[' '; s(:,1); ' ']';
     if ~nargout && ~isempty(istr), s={' ', [ '   ' istr ' ='], s{:}}; end

     s(2,:)={ char(10) }; s{2,end}='';

     s=regexprep(cat(2,s{:}),'--[ :=]*--','');
  end

  if ~nargout, disp(s); clear s; return; end

end

